
using MCMCChains, LinearAlgebra
include("Likelihood.jl")
include("Support.jl")


function PMMH(bootstrapFilter::Function, Y::DataFrame, opts::Dict; verbose=false)
    
    # Extract options
    nChains = opts["nChains"]
    chunkSize = opts["chunkSize"]
    maxChunks = opts["maxChunks"]
    nParams = length(opts["paramPriors"])
    windin = opts["windin"]
    covMat = opts["proposalCovariance"]
    
    maxSamples = maxChunks * chunkSize
    
    # Check that we are actually running multiple chainss
    if nChains < 2
        error("These methods require a minimum of two PMMH chains to be run.")
    end
    
    # Pre-allocate output
    θall = zeros(maxSamples, nParams, nChains) # A nSamples x nParams x nChains array
    diag = DataFrame() # A diagnostics dataframe, note the different format to simplePMMH()
    
    (ii, ESS, Rhat) = (1, 0, Inf)
    while (ii <= opts["maxChunks"] && (Rhat > opts["maxRhat"] || ESS < opts["minESS"]))
        
        # Fetch inds to save to
        firstInd = (ii-1)*chunkSize + 1
        lastInd = ii*chunkSize
        inds = firstInd:lastInd
        
        # Run a single chunk of PMMH for each chain
        nAcceptAll = zeros(nChains)
        ProgBar = Progress(nChains*(chunkSize-1), dt=1, desc="Running chunk $ii", barlen=50, enabled=verbose)
        Threads.@threads for jj = 1:nChains
            initialValues = missing
            if ii > 1
                initialValues = θall[firstInd-1,:,jj]
            end
            (θall[inds,:,jj], nAcceptAll[jj]) = runSingleEpochPMMH(bootstrapFilter, covMat, Y, opts; ProgBar=ProgBar, chainNum=jj, initialValues=initialValues)
        end
        finish!(ProgBar)
        
        if lastInd > windin
            
            chains = Chains(θall[windin:lastInd,:,:])
            dfc = DataFrame(summarize(chains))
            
            # Fetch max Rhat and min ESS
            Rhat = maximum(dfc.rhat)
            ESS = minimum(dfc.ess_bulk)
            
            # Save diagnostics
            diag_tmp = DataFrame()
            diag_tmp.chunk = repeat([ii], nParams)
            diag_tmp.pAccept = repeat([sum(nAcceptAll)/(nChains * (chunkSize - 1))],  nParams)
            diag_tmp.Rhat = dfc.rhat
            diag_tmp.ESS = dfc.ess_bulk
            diag_tmp.mean = dfc.mean
            diag_tmp.std = dfc.std
            diag = vcat(diag, diag_tmp)
            
            if verbose
                println(diag_tmp)
            end
            
        end
        
        # Increment the number of chunks
        ii = ii + 1
        
    end
    
    lastInd = (ii-1)*chunkSize
    return(θall[windin:lastInd, :, :], diag)
    
end









function runSingleEpochPMMH(bootstrapFilter::Function, covEst, Y::DataFrame, opts::Dict; ProgBar=missing, chainNum=0, initialValues=missing)
    
    # Extract options that are frequently used
    chunkSize = opts["chunkSize"]
    nParams = length(opts["paramPriors"])
    
    # Pre-allocate outputs
    θ = zeros(chunkSize, nParams)
    loglik = -Inf
    nAccept = 0
    
    # Sample initial parameter values and compute initial log-likelihood (only if no previous values are included)
    if ismissing(initialValues)
        (θ[1,:], loglik) = sampleInitialValues(bootstrapFilter, Y, opts)
    else
        θ[1,:] = initialValues
        loglik = estimateLogLik(bootstrapFilter, initialValues, Y, opts; ignoreerror=true)
    end
    
    # Setup progress bar if one is not provided
    if ismissing(ProgBar)
        ProgBar = Progress(chunkSize-1, dt=1, desc="Running single chunk...", barlen=50, enabled=opts["showEpochProgress"])
    end
    
    # Run PMMH
    for ii = 2:chunkSize
        
        # # Construct proposal distributions
        # proposalDists = [Truncated(Normal(θ[ii-1,jj], propStdDev[jj]), opts["paramLimits"][jj][1], opts["paramLimits"][jj][2]) for jj in 1:nParams]
        
        # # Propose new parameters and estimate log-likelihood
        # θprop = [rand(proposalDists[jj]) for jj in 1:nParams]
            θprop = rand(MvNormal(θ[ii-1,:], ((2.38^2)/nParams) * covEst))
            
            
            loglikProp = estimateLogLik(bootstrapFilter, θprop, Y, opts; ignoreerror=true)
            
            # Calculate acceptance probability
            # proposalDistsBackward = [Truncated(Normal(θprop[jj], propStdDev[jj]), opts["paramLimits"][jj][1], opts["paramLimits"][jj][2]) for jj in 1:nParams]
            α = calculateAcceptanceProbability(loglikProp, loglik, θprop, θ[ii-1,:], covEst, opts)
            
            # Accept or reject proposal
            if rand() < α
                θ[ii,:] = θprop
                loglik = loglikProp
                nAccept += 1
            else
                θ[ii,:] = θ[ii-1,:]
            end
            
            next!(ProgBar)
            
        end
        
        return(θ, nAccept)
        
    end
    
    
    
    function sampleInitialValues(bootstrapFilter::Function, Y::DataFrame, opts::Dict)
        
        initalValuesSampled = false
        nAttempts = 0
        θ = missing
        loglik = missing
        
        while !initalValuesSampled
            nAttempts = nAttempts + 1
            θ = [rand(opts["initialParamSamplers"][ii]) for ii in 1:length(opts["initialParamSamplers"])]
            try
                loglik = estimateLogLik(bootstrapFilter, θ, Y, opts)
                initalValuesSampled = true
            catch e
                println("Initial value θ = $θ failed. Trying again...")
            end
            if nAttempts > 100
                error("Error sampling initial values. Either initial values could not be sampled, or there is a bug in the code.")
            end
        end
        
        return(θ, loglik)
        
    end
    
    
    function calculateAcceptanceProbability(loglikprop, loglikprev, θprop, θprev, covEst, opts)
        
        nParams = length(θprop)
        
        # Likelihood ratio
        likelihoodratio = exp(loglikprop - loglikprev)
        
        priorratio = prod([pdf(opts["paramPriors"][jj], θprop[jj]) / pdf(opts["paramPriors"][jj], θprev[jj]) for jj in 1:nParams])
        
        # proposalProbsForward = [pdf(proposalDists[jj], θprop[jj]) for jj in 1:nParams]
        # proposalProbsBackward = [pdf(proposalDistsBackward[jj], θprev[jj]) for jj in 1:nParams]
        # proposalratio = prod(proposalProbsForward) / prod(proposalProbsBackward)
        proposalratio = pdf(MvNormal(θprev, covEst), θprop) / pdf(MvNormal(θprop, covEst), θprev)
        
        acceptanceprob = min(1, likelihoodratio * priorratio * proposalratio)
        
        return(acceptanceprob)
        
    end
    
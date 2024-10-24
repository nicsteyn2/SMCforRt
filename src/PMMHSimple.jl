
using MCMCChains, CSV
include("Likelihood.jl")

# --------------------------------------------------------------------------
# ------------------------- Core PMMH functions ----------------------------
# --------------------------------------------------------------------------

function multipleSimplePMMH(bootstrapFilter::Function, Y::DataFrame, opts::Dict; showProgress=true)

    # Check that we are actually running multiple
    if opts["nChains"] == 1
        error("Please set nChains > 1 to use runMultipleSimplePMMH(), otherwise use simplePMMH().")
    end

    # Pre-allocate output
    θall = zeros(opts["nPMMHSamples"], length(opts["paramPriors"]), opts["nChains"])
    diagnosticsAll = Array{Dict}(undef, opts["nChains"])

    ProgBar = Progress((opts["nPMMHSamples"]-1)*opts["nChains"], dt=1, desc="Running multiple PMMH...", barlen=50, enabled=showProgress)
    Threads.@threads for ii = 1:opts["nChains"]
        (θall[:,:,ii], diagnosticsAll[ii]) = simplePMMH(bootstrapFilter, Y, opts, ProgBar=ProgBar, chainNum=ii)
    end

    return(θall, diagnosticsAll)
    
end


function simplePMMH(bootstrapFilter::Function, Y::DataFrame, opts::Dict; ProgBar=missing, chainNum=0, showProgress=true)

    # Extract options that are frequently used
    nSamples = opts["nPMMHSamples"]
    nParams = length(opts["paramPriors"])

    # Pre-allocate outputs
    θ = zeros(nSamples, nParams)
    logliks = zeros(nSamples)
    nAccept = 0

    # Sample initial parameter values and compute initial log-likelihood
    θ[1,:] = [rand(opts["initialParamSamplers"][ii]) for ii in 1:nParams]
    logliks[1] = estimateLogLik(bootstrapFilter, θ[1,:], Y, opts)

    # Create progress bar if one hasn't been provided
    if ismissing(ProgBar)
        ProgBar = Progress(nSamples-1, dt=1, desc="Running PMMH...", barlen=50, enabled=showProgress)
    end

    # Run PMMH
    for ii = 2:nSamples
        
        # Propose new parameters, estimate log-likelihood, and store the sample
        θprop = [rand(opts["proposalDists"][jj](θ[ii-1,jj])) for jj in 1:nParams]
        loglikProp = NaN
        try
            loglikProp = estimateLogLik(bootstrapFilter, θprop, Y, opts)
        catch e
            loglikProp = -Inf
        end


        # Compute acceptance probability
        likelihoodratio = exp(loglikProp - logliks[ii-1])

        priorratio = prod([pdf(opts["paramPriors"][jj], θprop[jj]) / pdf(opts["paramPriors"][jj], θ[ii-1,jj]) for jj in 1:nParams])

        proposalProbsForward = [pdf(opts["proposalDists"][jj](θ[ii-1,jj]), θprop[jj]) for jj in 1:nParams]
        proposalProbsBackward = [pdf(opts["proposalDists"][jj](θprop[jj]), θ[ii-1,jj]) for jj in 1:nParams]
        proposalratio = prod(proposalProbsForward) / prod(proposalProbsBackward)

        acceptanceprob = min(1, likelihoodratio * priorratio * proposalratio)

        # Accept or reject proposal
        if rand() < acceptanceprob
            θ[ii,:] = θprop
            logliks[ii] = loglikProp
            nAccept += 1
        else
            θ[ii,:] = θ[ii-1,:]
            logliks[ii] = logliks[ii-1]
        end

        # Update progress bar
        next!(ProgBar)

    end

    finish!(ProgBar)

    # Store diagnostics and return
    diagnostics = Dict()
    diagnostics["logliks"] = logliks
    diagnostics["pAccept"] = nAccept/(nSamples-1)

    return(θ, diagnostics)

end






# --------------------------------------------------------------------------
# ---------------------------- Helper functions ----------------------------
# --------------------------------------------------------------------------


# This calculates the theoretical optimal standard deviation of the proposal distribution and compares
# it to the actual standard deviation used in the PMMH algorithm.
function getProposalStdDevTable(chains, opts)

    nParams = sum(opts["proposalStdDev"] .> 0)
    posteriorStdDev = DataFrame(summarystats(chains)).std
    out = DataFrame(propStd=proposalStdDevs, theoreticalStd=(2.38*posteriorStdDev)/nParams)
    out.relative = out.theoreticalStd ./ out.propStd
    return(out)

end


function savePMMHSummary(θ, didAccept, opts; print=true, windin=100, paramNames=false)

    fname = opts["filename"] * "summary.txt"

    open(fname, "w") do sink

        # Print header
        println(sink, "---------------------------------------------------------------")
        println(sink, "------------------- Summary of PMMH results -------------------")
        println(sink, "---------------------------------------------------------------")
        println(sink, "")

        # Print chains
        println(sink, "Removing first $(windin) samples from each chain")
        println(sink, "")
        if paramNames != false
            if (length(paramNames) >= 2) && (paramNames[2] == "ρ")
                θtemp = deepcopy(θ)
                θtemp[:,2,:] = θtemp[:,2,:]*1e4
                show(sink, MIME"text/plain"(), Chains(θtemp[(windin+1):end,:,:], paramNames)) 
            elseif (length(paramNames) >= 6) && (paramNames[6] == "ρ")
                θtemp = deepcopy(θ)
                θtemp[:,6,:] = θtemp[:,6,:]*1e4
                show(sink, MIME"text/plain"(), Chains(θtemp[(windin+1):end,:,:], paramNames))
            else
                show(sink, MIME"text/plain"(), Chains(θ[(windin+1):end,:,:], paramNames))
            end
        else
            show(sink, MIME"text/plain"(), Chains(θ[(windin+1):end,:,:]))
        end
        println(sink, "")

        # Print acceptance rates
        println(sink, "")
        println(sink, "---------------------------------------------------------------")
        println(sink, "Overall acceptance probability = ", mean(didAccept), ". By chain:")
        println(sink, mean(didAccept, dims=1))
        println(sink, "---------------------------------------------------------------")
        println(sink, "")

        # Print proposal std dev table
        propStdDevTable = getProposalStdDevTable(Chains(θ[(windin+1):end,:,:]), opts)
        if paramNames != false
            propStdDevTable = hcat(DataFrame(param=paramNames), propStdDevTable)
        end
        # Print propStdDevTable but don't print row numbers
        println(sink, "---------------------------------------------------------------")
        println(sink, "Theoretical vs actual proposal standard deviation")
        println(sink, propStdDevTable)
        println(sink, "---------------------------------------------------------------")
        println(sink, "")

    end

    if print
        println(String(read(fname)))
    end

end
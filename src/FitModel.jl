
using DataFrames

include("Support.jl")
include("PMMH.jl")
include("MarginalPosterior.jl")
include("Likelihood.jl")

# This performs a general start-to-finish model fitting procedure

function fitModel(bootstrapFilter::Function, Y::DataFrame, optsIn::Dict; skipResamplingPMMH=false, predictiveValues=missing, verbose=false, checkStdDevLogLik=true)  #TOOD: Include predcitve values fucnctionality
    
    # Ensure we do not edit the opts dictionary
    opts = deepcopy(optsIn)
    
    # Extract options
    paramNames = opts["paramNames"] # PMMH parameter names (useful for storing tidy output)
    stateNames = opts["stateNames"] # Hidden-state names
    nParams = length(opts["paramPriors"])
    θtest = missing
    if haskey(opts, "initialParamSamplers")
        θtest = [mean(opts["initialParamSamplers"][ii]) for ii in 1:nParams]
    else
        θtest = [mean(opts["paramPriors"][ii]) for ii in 1:nParams]
    end
    
    # Start the timer
    S = time()
    
    # If the model is Markovian, we can skip resampling when estimating the likelihood
    L = opts["L"] # Store this to use when finding the marginal posterior distribution
    if skipResamplingPMMH
        opts["L"] = 1
    end
    
    # Check the standard deviation of the log-likelihood
    if checkStdDevLogLik
        (loglikstddev, _) = estimateStdDevLogLik(100, bootstrapFilter, θtest, Y, opts; showProgress=verbose)
        if verbose
            println("Estimated standard deviation of the log-likelihood: ", loglikstddev)
        end
        if loglikstddev > 4
            error("The standard deviation of the log-likelihood is too high (>4). Consider increasing the number of particles used, or changing the model.")
        elseif loglikstddev > 2
            @warn("The standard deviation of the log-likelihood is high (>2). Consider increasing the number of particles used, or changing the model.")
        end
    end
    
    # Run PMMH
    (θ, diag) = PMMH(bootstrapFilter, Y, opts; verbose=verbose)
    chains = Chains(θ, opts["paramNames"])
    
    # Reset resampling window
    opts["L"] = L
    
    # Store parameter values
    df_params = DataFrame()
    df_params.param = paramNames
    df_params.m = [mean(θ[:,ii,:]) for ii in 1:size(θ)[2]]
    df_params.l = [quantile(vec(θ[:,ii,:]), 0.025) for ii in 1:size(θ)[2]]
    df_params.u = [quantile(vec(θ[:,ii,:]), 0.975) for ii in 1:size(θ)[2]]
    df_params.rhat = DataFrame(summarize(chains)).rhat
    df_params.ess = DataFrame(summarize(chains)).ess
    
    # Fetch samples from the marginal posterior distribution
    X = marginalPosterior(bootstrapFilter, θ, Y, opts)
    
    # Extract dimensions
    nStates = missing
    if length(size(X)) == 2
        nStates = 1
    else
        nStates = size(X)[3]
    end
    
    # Process marginal posterior samples
    df_states = DataFrame()
    if nStates == 1
        (m, _, l, u) = processResults(X)
        tmp = DataFrame(t=Y.t, mean=m, lower=l, upper=u, variable=stateNames[1])
        df_states = vcat(df_states, tmp)
    else
        for ii = 1:nStates
            (m, _, l, u) = processResults(X[:,:,ii])
            tmp = DataFrame(t=Y.t, mean=m, lower=l, upper=u, variable=stateNames[ii])
            df_states = vcat(df_states, tmp)
        end
    end
    
    # If Y contains a date vector, include it in our outputs
    if "date" in names(Y)
        df_states = leftjoin(df_states, Y[:,[:t, :date]], on=:t)
    end
    
    # End the timer
    E = time()
    if verbose
        println("Model finished. Time elapsed: ", round(E-S, sigdigits=3), " seconds")
    end
    
    return(df_states, df_params, θ, diag)
    
end





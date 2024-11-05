
using ProgressMeter

function marginalPosterior(bootstrapFilter::Function, θ, Y::DataFrame, opts::Dict; showProgress=false) #TOOD: Include predcitve values fucnctionality
    
    # Make a copy of the options dictionary so we can edit
    optsIn = deepcopy(opts)
    optsIn["N"] = optsIn["posteriorNumberOfParticles"]
    
    # Extract frequently used options
    N = optsIn["N"]
    nPosteriorParamSamples = optsIn["posteriorParamSamples"]
    
    # Calculate derived quantities
    Nout = nPosteriorParamSamples * optsIn["posteriorNumberOfParticles"]
    
    # If θ is passed in directly from PMMH, resize it
    if length(size(θ)) == 3
        θ = resizeParams(θ)
    elseif length(size(θ)) != 2
        error("Invalid parameter matrix dimensions.")
    end
    
    # Sample parameter values
    θsample = θ[sample(1:size(θ)[1], nPosteriorParamSamples, replace=false),:]
    
    # Run the filter once to determine size of output
    (X0, _) = bootstrapFilter(θsample[1,:], Y, optsIn)
    
    # Extract dimensions
    nStates = missing
    if length(size(X0)) == 2
        T = size(X0)[2]
        nStates = 1
    else
        T = size(X0)[2]
        nStates = size(X0)[3]
    end
    
    # Pre-allocate output
    X = zeros(Nout, T, nStates)
    
    # Store first sample
    if nStates == 1
        X[1:N,:,1] = X0
    else
        for ii = 1:nStates
            X[1:N,:,ii] = X0[:,:,ii]
        end
    end
    
    # Run the model at the remaining parameters
    ProgBar = Progress(nPosteriorParamSamples-1, dt=1, barlen=50, desc="Sampling from marginal posterior...", enabled=showProgress)
    Threads.@threads for ii = 2:nPosteriorParamSamples
        (X0, _) = bootstrapFilter(θsample[ii,:], Y, optsIn)
        inds = ((ii-1)*N + 1):(ii*N)
        if nStates == 1
            X[inds,:,1] = X0
        else
            for jj = 1:nStates
                X[inds,:,jj] = X0[:,:,jj]
            end
        end
        next!(ProgBar)
    end

    # Resize if we are working with a single state
    if nStates == 1
        X = X[:,:,1]
    end
    
    return(X)
    
end



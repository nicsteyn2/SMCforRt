
using ProgressMeter

function estimateLogLik(bootstrapFilter::Function, θ, Y::DataFrame, opts::Dict; ignoreerror=false)

    try
        (_, W) = bootstrapFilter(θ, Y, opts)
        loglik = sum(log.(mean(W, dims=1)[2:end]))
        return(loglik)
    catch e
        if ignoreerror
            return(-Inf)
        else
            rethrow(e)
        end
    end

end
    



function estimateStdDevLogLik(nTrials::Int, bootstrapFilter::Function, θ, Y::DataFrame, opts::Dict; showProgress=true, ignoreerror=false)

    logliks = zeros(nTrials)
    ProgBar = Progress(nTrials, dt=1, desc="Estimating log-likelihood multiple times...", barlen=50, enabled=showProgress)
    Threads.@threads for ii in 1:nTrials
        logliks[ii] = estimateLogLik(bootstrapFilter, θ, Y, opts; ignoreerror=ignoreerror)
        next!(ProgBar)
    end
    
    return(std(logliks), logliks)

end
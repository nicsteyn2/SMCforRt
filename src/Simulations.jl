
using Distributions, DataFrames

include("Support.jl") # Contains the default generation time distribution

function simulateSimpleEpidemic(; σ=0.1, g=missing, R1=2, I1=10, T=100, minInfections=100, maxInfections=100000, minRt = 0.2, maxRt = 5)
    
    attemptsRemaining = 100
    isFinished = false
    
    # If we do not specify a generation time distribution, we use the default
    if ismissing(g)
        g = defaultGenTimeDist(T)
    end
    
    # Pre-allocate output
    logRt = zeros(T)
    logRt[1] = log(R1)
    It = Int.(zeros(T))
    It[1] = I1
    
    while !isFinished && (attemptsRemaining > 0)
        
        It = Int.(zeros(T))
        It[1] = I1
        
        for tt = 2:T
            logRt[tt] = rand(Truncated(Normal(logRt[tt-1], σ), log(minRt), log(maxRt))) # Sample from a truncated normal distribution
            Λt = sum(It[tt-1:-1:1] .* g[1:tt-1]) # Calculate the force of infection
            It[tt] = rand(Poisson(exp(logRt[tt]) * Λt)) # Sample the number of new infections
            if sum(It[1:tt]) > maxInfections
                break
            end
        end
        
        if (minimum(It) >= minInfections) && (maximum(It) <= maxInfections)
            isFinished = true
        end

        attemptsRemaining -= 1
        
    end
        
    if attemptsRemaining == 0
        error("Failed to simulate the epidemic in the allowed number of attempts. Check min and max infection requirements.")
    else
        return(DataFrame(t = 1:T, It = It, Rt = exp.(logRt)))
    end
    
end





# Simulate observation processes

function simulateBinomialReporting(It, ρ)

    Ct = rand.(Binomial.(It, ρ))
    return(Ct)

end



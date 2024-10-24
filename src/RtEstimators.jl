
function CalculateInfectionPressure(w::Vector, C::Vector; normalised=true)
    
    Λ = missing
    if normalised
        Λ = [sum(C[(tt-1):-1:1] .* w[1:(tt-1)])/sum(w[1:max(tt-1, 1)]) for tt = 1:length(C)]
    else
        Λ = [sum(C[(tt-1):-1:1] .* w[1:(tt-1)]) for tt = 1:length(C)]
    end
    return(Λ)
    
end




# --------------------------------------------------
# ------------------- EpiEstim ---------------------
# --------------------------------------------------

function EpiEstimPosteriorParams(k::Int, w::Vector, C::Vector; a0=1, b0=1/5)
    
    # Calculate infection pressure
    Λ = CalculateInfectionPressure(w, C)
    
    # Preallocate output
    a = zeros(length(C))
    b = zeros(length(C))
    
    # Run for each step (careful about parameters for tt<k)
    if k == 0
        a[1:length(C)] .= a0
        b[1:length(C)] .= b0
    else
        a[1] = a0
        b[1] = b0 # At t = 1 there is no past data to inform b
        for tt = 2:length(C)
            a[tt] = a0 + sum(C[max(tt-k+1,1):tt])
            b[tt] = b0 + sum(Λ[max(tt-k+1,1):tt])
        end
    end
    
    return(a, b)
    
end

function EpiEstim(k::Int, w::Vector, C::Vector; a0=1, b0=1/5)
    
    # Calculate infection pressure
    Λ = CalculateInfectionPressure(w::Vector, C::Vector)
    
    # Fetch posterior parameters and posterior dist
    (a, b) = EpiEstimPosteriorParams(k, w, C; a0=a0, b0=b0)
    pRgivenk = Gamma.(a, 1 ./ b)
    
    return(pRgivenk)
    
end






function CRPS(X::Vector, y)

    N = length(X)
    term1 = mean(abs.(X .- y))
    term2 = 0.5 * (1/N^2) * sum( abs.(X[i] - X[j]) for i in 1:N, j in 1:N)
    return(term1 - term2)

end


function CRPS(X::Matrix, y::Vector)

    nCols = size(X, 2)
    CRPS_col = zeros(nCols)
    for ii = 1:nCols
        CRPS_col[ii] = CRPS(X[:,ii], y[ii])
    end
    return(mean(CRPS_col))

end
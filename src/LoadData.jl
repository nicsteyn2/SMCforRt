
using DataFrames, CSV, Dates, Distributions


# Loads various data-sources.
# Returns a dataframe typically containing :t (time) and :Ct (cases), but this may vary.
function loadData(fname)

    if fname == "simulated_simple"

        return CSV.read("../data/simulated_simple.csv", DataFrame)

    elseif fname == "NZCOVID" # Data obtained from https://github.com/minhealthnz/nz-covid-data

        nzdata = CSV.read("../data/nzcovid_moh.csv", DataFrame)
        sort!(nzdata, :date)
        nzdata.Ct = nzdata.border + nzdata.local
        nzdata.t = Dates.value.(nzdata.date .- minimum(nzdata.date)) .+ 1
        return nzdata[1:100,:]
        
    elseif fname == "NZCOVID_1APR2024"

        nzdata = CSV.read("../data/nzcovid_moh.csv", DataFrame)
        sort!(nzdata, :date)
        nzdata.Ct = nzdata.border + nzdata.local
        nzdata.t = Dates.value.(nzdata.date .- minimum(nzdata.date)) .+ 1
        return nzdata[1497:end,:]
        
    end

end
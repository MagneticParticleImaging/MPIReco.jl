using HTTP

if !isdir("dataMP01")
  HTTP.open("GET", "http://media.tuhh.de/ibi/mpireco/data.zip") do http
    open("data.zip", "w") do file
        write(file, http)
    end
  end
  run(`unzip data.zip`)
end

filenameSM = "systemMatrix.mdf"
filenameMeas = "measurement.mdf"

if !isfile(filenameSM)
  HTTP.open("GET", "http://media.tuhh.de/ibi/mdfv2/systemMatrix_V2.mdf") do http
    open(filenameSM, "w") do file
        write(file, http)
    end
  end
end
if !isfile(filenameMeas)
  HTTP.open("GET", "http://media.tuhh.de/ibi/mdfv2/measurement_V2.mdf") do http
    open(filenameMeas, "w") do file
        write(file, http)
    end
  end
end


mkpath("./img/")

include("MultiPatch.jl")
include("MultiGradient.jl")
include("Regular.jl")
include("Sparse.jl")
include("SMCenter.jl")

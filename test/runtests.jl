import ImageMagick
using HTTP
using Test
using FileIO

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

mkpath("./motionComp/")
for f in ["SF1Small.mdf", "SF2Small.mdf", "SF3Small.mdf", "SF4Small.mdf",
          "measBG.mdf", "measFast.mdf"] #, "measSlow.mdf"]
  filename = joinpath("./motionComp/", f)
  if !isfile(filename)
    HTTP.open("GET", joinpath("http://media.tuhh.de/ibi/motionComp/",f)) do http
      open(filename, "w") do file
          write(file, http)
      end
    end
  end

end

mkpath("./img/")

function exportImage(filename, I::AbstractMatrix)
  Iabs = abs.(I)
  Icolored = colorview(Gray, Iabs./maximum(Iabs))
  save(filename, Icolored )
end

include("Reconstruction.jl")
include("MultiPatch.jl")
include("MultiGradient.jl")
include("Sparse.jl")
include("SMCenter.jl")
include("LoadSaveMDF.jl")
include("SMRecovery.jl")
include("MotionCompensation.jl")

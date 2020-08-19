import ImageMagick
using HTTP
using Test
using FileIO

if !isdir("data")
  @info "download data.zip"
  HTTP.open("GET", "http://media.tuhh.de/ibi/MPIReco/data.zip") do http
    open("data.zip", "w") do file
        write(file, http)
    end
  end
  @info "extracting data.zip"
  run(`unzip -oq data.zip`)
  rm("data.zip")
end

mkpath("./img/")

function exportImage(filename, I::AbstractMatrix)
  Iabs = abs.(I)
  Icolored = colorview(Gray, Iabs./maximum(Iabs))
  save(filename, Icolored )
end

include("Cartesian.jl")
include("MotionCompensation.jl")
include("Reconstruction.jl")
include("MultiPatch.jl")
include("MultiGradient.jl")
include("Sparse.jl")
include("SMCenter.jl")
include("LoadSaveMDF.jl")
include("SMRecovery.jl")

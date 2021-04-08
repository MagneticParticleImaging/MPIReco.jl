import ImageMagick
using HTTP
using Test
using FileIO
using Scratch
using LazyArtifacts

const datadir = joinpath(artifact"data", "data")
const imgdir  = @get_scratch!("img")

@info "If you want to check the output of the tests, please head to $imgdir"

function exportImage(filename, I::AbstractMatrix)
  Iabs = abs.(I)
  Icolored = colorview(Gray, Iabs./maximum(Iabs))
  save(filename, Icolored )
end

include("LoadSaveMDF.jl")
include("Reconstruction.jl")
include("Cartesian.jl")
include("MotionCompensation.jl")
include("MultiPatch.jl")
include("MultiGradient.jl")
include("Sparse.jl")
include("SMCenter.jl")
include("SMRecovery.jl")

using HTTP
using Test
using FileIO
using Scratch
using LazyArtifacts

const datadir = joinpath(artifact"data")
@info "The test data is located at $datadir."

const imgdir  = @get_scratch!("img")
@info "If you want to check the output of the tests, please head to $imgdir."

function exportImage(filename, I::AbstractMatrix)
  Iabs = abs.(I)
  Icolored = colorview(Gray, Iabs./maximum(Iabs))
  save(filename, Icolored )
end

#include("LoadSaveMDF.jl")
#include("Reconstruction.jl")
#include("Cartesian.jl")
#if !Sys.iswindows()
#  include("MotionCompensation.jl")
#end
#include("MultiPatch.jl")
#include("MultiGradient.jl")
#include("Sparse.jl")
#include("SMCenter.jl")
#include("SMRecovery.jl")
include("SMExtrapolation.jl")
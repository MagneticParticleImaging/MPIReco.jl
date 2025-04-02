using MPIReco
using AbstractImageReconstruction
using JLArrays
using Test
using Scratch
using ImageMagick
using ImageQualityIndexes
using LazyArtifacts
using Unitful
using Dagger
using DaggerImageReconstruction

const datadir = joinpath(artifact"data")
@info "The test data is located at $datadir."

const imgdir  = @get_scratch!("img")
@info "If you want to check the output of the tests, please head to $imgdir."

function compareImg(filename, cmp::Function = compareSSIM, kwargs...)
  @info "Comparing image $filename"
  imExpected = ImageMagick.load(joinpath("correct", filename));
  imGiven = ImageMagick.load(joinpath(imgdir, filename));
  return cmp(imExpected, imGiven, kwargs...)
end

const SSIM_LIMIT = 0.97
function compareSSIM(expected, given; limit::Float64=SSIM_LIMIT, kwargs...)
  ssim = assess_ssim(given, expected)
  @info "Image Quality (SSIM): $ssim"
  return ssim >= limit
end

function exportImage(filename, I::AbstractMatrix)
  Iabs = abs.(I)
  Icolored = colorview(Gray, Iabs./maximum(Iabs))
  ImageMagick.save(filename, Icolored )
end

areTypesDefined = @isdefined arrayTypes
arrayTypes = areTypesDefined ? arrayTypes : [JLArray]

@testset "MPIReco" begin
  if !areTypesDefined
    include("Caching.jl")
    #include("LoadSaveMDF.jl")
    include("Reconstruction.jl") # FussedLasso causes segfault atm
    include("Solvers.jl")
    include("Cartesian.jl")
    if !Sys.iswindows()
     include("MotionCompensation.jl")
    end
    include("MultiPatch.jl")
    include("MultiGradient.jl")
    include("Sparse.jl")
    include("SMExtrapolation.jl")
    include("ReconstructionGPU.jl")
    #include("SMCenter.jl")
    #include("SMRecovery.jl") # will be moved to SMTools
  else
    include("ReconstructionGPU.jl")
  end
end

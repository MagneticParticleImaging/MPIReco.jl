import ImageMagick
using HTTP
using Test
using FileIO

using Scratch
using Artifacts

# TODO: Don't use ZIP, instead upload a tar.gz to somewhere (like Github)
#       https://github.com/JuliaPackaging/Scratch.jl#can-i-use-a-scratch-space-as-a-temporary-workspace-then-turn-it-into-an-artifact

# To create a new Artifact.toml use https://github.com/simeonschaub/ArtifactUtils.jl

scratch = @get_scratch!("data")
if isempty(readdir(scratch))

  @info "download data"
  datapath = joinpath(artifact"data", "data.zip")

  @info "extracting data"
  cd(scratch) do
    run(`unzip -oq $datapath`)
  end
end

const datadir = joinpath(scratch, "data")
const imgdir  = @get_scratch!("img")

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

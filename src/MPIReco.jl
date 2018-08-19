module MPIReco

using Pkg

if !haskey(Pkg.installed(),"LinearSolver")
  println("Installing LinearSolver...")
  Pkg.clone("https://github.com/tknopp/LinearSolver.jl.git")
end

if !haskey(Pkg.installed(),"MPIFiles")
  println("Installing MPIFiles...")
  Pkg.clone("https://github.com/MagneticParticleImaging/MPIFiles.jl.git")
end

using Reexport
@reexport using MPIFiles
@reexport using LinearSolver
@reexport using Images
using AxisArrays
const axes = Base.axes
using Compat
using ProgressMeter
using LinearAlgebra
using Statistics
using Random
using SparseArrays
using Dates

include("Utils.jl")
include("RecoParameters.jl")
include("SystemMatrixCenter.jl")
include("SystemMatrix.jl")
include("Weighting.jl")
include("Reconstruction.jl")
include("MultiPatch.jl")


end # module

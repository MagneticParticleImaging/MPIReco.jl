module MPIReco

if !isdir(Pkg.dir("LinearSolver"))
  println("Installing LinearSolver...")
  Pkg.clone("https://github.com/tknopp/LinearSolver.jl.git")
end

if !isdir(Pkg.dir("MPIFiles"))
  println("Installing MPIFiles...")
  Pkg.clone("https://github.com/MagneticParticleImaging/MPIFiles.jl.git")
end

using Reexport
@reexport using MPIFiles
@reexport using LinearSolver
@reexport using Images
using Compat
using ProgressMeter

include("Utils.jl")
include("SystemMatrix.jl")
include("Weighting.jl")
include("Reconstruction.jl")
include("MultiPatch.jl")


end # module

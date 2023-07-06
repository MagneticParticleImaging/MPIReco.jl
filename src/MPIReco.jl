module MPIReco
  using Reexport
  @reexport using RegularizedLeastSquares
  @reexport using ImageUtils
  @reexport using MPIFiles
  const shape = MPIFiles.shape
  using RecoUtils
  @reexport using DSP
  using ProgressMeter
  using LinearAlgebra
  using SparseArrays
  using Distributed
  using DistributedArrays
  # using TensorDecompositions
  using IniFile
  import LinearAlgebra: ldiv!, \
  # TODO sort out import for Base and RecoUtils to avoid boiler plate
  import Base: put!, take!

  # ML Stuff
  import Distributions
  using Flux
  using CUDA
  using BSON

  # SMextrapolation
  import SparseArrays.spdiagm

  include("AlgorithmInterface.jl")
  include("Background.jl")
  include("PreProcessing/PreProcessing.jl")
  include("LeastSquares.jl")
  include("Utils.jl")
  include("MultiContrast.jl")
  include("RecoParameters.jl")
  include("SystemMatrix/SystemMatrix.jl")
  include("Weighting.jl")
  include("TemporalRegularization/TemporalRegularization.jl")
  include("Reconstruction.jl")
  include("MultiPatch.jl")
  include("MotionCompensation/MotionCompensation.jl")
  include("MachineLearning/MachineLearning.jl")
  include("Algorithms/Algorithms.jl")
  include("Serialisation.jl")
end # module

module MPIReco
  using RegularizedLeastSquares
  @reexport using RegularizedLeastSquares
  @reexport using ImageUtils
  @reexport using MPIFiles
  const shape = MPIFiles.shape
  @reexport using DSP

  using LinearAlgebra
  using SparseArrays
  using Distributed
  using DistributedArrays
  # using TensorDecompositions
  using IniFile
  import LinearAlgebra: ldiv!, \

  # ML Stuff
  import Distributions
  using Flux
  using CUDA
  using BSON


  include("Utils.jl")
  include("MultiContrast.jl")
  include("RecoParameters.jl")
  include("SystemMatrix/SystemMatrix.jl")
  include("Weighting.jl")
  include("Background.jl")
  include("TemporalRegularization/TemporalRegularization.jl")
  include("Reconstruction.jl")
  include("MultiPatch.jl")
  include("MotionCompensation/MotionCompensation.jl")
  include("MachineLearning/MachineLearning.jl")
end # module

module MPIReco
  using Reexport
  @reexport using RegularizedLeastSquares
  using RegularizedLeastSquares.LinearOperators
  @reexport using ImageUtils
  @reexport using MPIFiles
  const shape = MPIFiles.shape
  @reexport using AbstractImageReconstruction
  using LRUCache
  using Adapt
  @reexport using DSP
  using ProgressMeter
  using Unitful
  using Statistics

  using LinearAlgebra
  using SparseArrays
  using Distributed
  using DistributedArrays
  # using TensorDecompositions
  using IniFile
  import LinearAlgebra: ldiv!, \, mul!
  # TODO sort out import for Base and AbstractImageReconstruction to avoid boiler plate
  import Base: put!, take!
  import AbstractImageReconstruction: process, parameter, reconstruct
  using FFTW
  using LinearOperatorCollection
  using RelocatableFolders


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
  include("Algorithms/Algorithms.jl")
  include("Serialisation.jl")
end # module

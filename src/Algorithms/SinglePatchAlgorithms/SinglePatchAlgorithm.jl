Base.@kwdef struct SinglePatchReconstructionParameter{L<:AbstractSystemMatrixLoadingParameter, SL<:AbstractLinearSolver,
   SP<:AbstractSolverParameters{SL}, R<:AbstractRegularization, W<:AbstractWeightingParameters} <: AbstractSinglePatchReconstructionParameters
  # File
  sf::MPIFile
  sfLoad::Union{L, ProcessResultCache{L}}
  # Solver
  solverParams::SP
  reg::Vector{R} = AbstractRegularization[]
  weightingParams::W = NoWeightingParameters()
end

Base.@kwdef mutable struct SinglePatchReconstructionAlgorithm{P} <: AbstractSinglePatchReconstructionAlgorithm where {P<:AbstractSinglePatchAlgorithmParameters}
  params::P
  # Could also do reconstruction progress meter here
  sf::Union{MPIFile, Vector{MPIFile}}
  S::AbstractArray
  grid::RegularGridPositions
  freqs::Vector{CartesianIndex{2}}
  output::Channel{Any}
end

function SinglePatchReconstruction(params::SinglePatchParameters{<:AbstractMPIPreProcessingParameters, R, PT}) where {R<:AbstractSinglePatchReconstructionParameters, PT <:AbstractMPIPostProcessingParameters}
  return SinglePatchReconstructionAlgorithm(params)
end
function SinglePatchReconstructionAlgorithm(params::SinglePatchParameters{<:AbstractMPIPreProcessingParameters, R, PT}) where {R<:AbstractSinglePatchReconstructionParameters, PT <:AbstractMPIPostProcessingParameters}
  freqs, S, grid = prepareSystemMatrix(params.reco)
  return SinglePatchReconstructionAlgorithm(params, params.reco.sf, S, grid, freqs, Channel{Any}(Inf))
end
recoAlgorithmTypes(::Type{SinglePatchReconstruction}) = SystemMatrixBasedAlgorithm()
AbstractImageReconstruction.parameter(algo::SinglePatchReconstructionAlgorithm) = algo.params

function prepareSystemMatrix(reco::SinglePatchReconstructionParameter{L,S}) where {L<:AbstractSystemMatrixLoadingParameter, S<:AbstractLinearSolver}
  freqs, sf, grid = process(AbstractMPIRecoAlgorithm, reco.sfLoad, reco.sf)
  sf, grid = prepareSF(S, sf, grid) 
  return freqs, sf, grid
end


AbstractImageReconstruction.take!(algo::SinglePatchReconstructionAlgorithm) = Base.take!(algo.output)

function process(algo::SinglePatchReconstructionAlgorithm, params::Union{A, ProcessResultCache{<:A}}, f::MPIFile, args...) where A <: AbstractMPIPreProcessingParameters
  result = process(typeof(algo), params, f, args...)
  if eltype(algo.S) != eltype(result)
    @warn "System matrix and measurement have different element data type. Mapping measurment data to system matrix element type."
    result = map(eltype(algo.S),result)
  end
  result = copyto!(similar(algo.S, size(result)...), result)
  return result
end


function process(algo::SinglePatchReconstructionAlgorithm, params::SinglePatchReconstructionParameter, u)
  weights = process(algo, params.weightingParams, u)

  B = getLinearOperator(algo, params)

  solver = LeastSquaresParameters(op = B, S = algo.S, reg = params.reg, solverParams = params.solverParams, weights = weights)

  result = process(algo, solver, u)

  return gridresult(result, algo.grid, algo.sf)
end

process(algo::SinglePatchReconstructionAlgorithm, params::Union{ChannelWeightingParameters, WhiteningWeightingParameters}, u) = map(real(eltype(algo.S)), process(typeof(algo), params, algo.freqs))

function getLinearOperator(algo::SinglePatchReconstructionAlgorithm, params::SinglePatchReconstructionParameter{<:DenseSystemMatixLoadingParameter, S}) where {S}
  return nothing
end

function getLinearOperator(algo::SinglePatchReconstructionAlgorithm, params::SinglePatchReconstructionParameter{<:SparseSystemMatrixLoadingParameter, S}) where {S}
  return createLinearOperator(params.sfLoad.sparseTrafo, eltype(algo.S); shape=tuple(shape(algo.grid)...))
end
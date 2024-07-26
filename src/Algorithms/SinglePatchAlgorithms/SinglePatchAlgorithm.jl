Base.@kwdef struct SinglePatchReconstructionParameter{L<:AbstractSystemMatrixLoadingParameter, SL<:AbstractLinearSolver,
   matT <: AbstractArray, SP<:AbstractSolverParameters{SL}, R<:AbstractRegularization, W<:AbstractWeightingParameters} <: AbstractSinglePatchReconstructionParameters
  # File
  sf::MPIFile
  sfLoad::Union{L, ProcessResultCache{L}}
  arrayType::Type{matT} = Array
  # Solver
  solverParams::SP
  reg::Vector{R} = AbstractRegularization[]
  weightingParams::Union{W, ProcessResultCache{W}} = NoWeightingParameters()
end

Base.@kwdef mutable struct SinglePatchReconstructionAlgorithm{P, matT <: AbstractArray} <: AbstractSinglePatchReconstructionAlgorithm where {P<:AbstractSinglePatchAlgorithmParameters}
  params::P
  # Could also do reconstruction progress meter here
  sf::Union{MPIFile, Vector{MPIFile}}
  S::AbstractArray
  arrayType::Type{matT}
  grid::RegularGridPositions
  freqs::Vector{CartesianIndex{2}}
  output::Channel{Any}
end

function SinglePatchReconstruction(params::SinglePatchParameters{<:AbstractMPIPreProcessingParameters, R, PT}) where {R<:AbstractSinglePatchReconstructionParameters, PT <:AbstractMPIPostProcessingParameters}
  return SinglePatchReconstructionAlgorithm(params)
end
function SinglePatchReconstructionAlgorithm(params::SinglePatchParameters{<:AbstractMPIPreProcessingParameters, R, PT}) where {R<:AbstractSinglePatchReconstructionParameters, PT <:AbstractMPIPostProcessingParameters}
  freqs, S, grid, arrayType = prepareSystemMatrix(params.reco)
  return SinglePatchReconstructionAlgorithm(params, params.reco.sf, S, arrayType, grid, freqs, Channel{Any}(Inf))
end
recoAlgorithmTypes(::Type{SinglePatchReconstruction}) = SystemMatrixBasedAlgorithm()
AbstractImageReconstruction.parameter(algo::SinglePatchReconstructionAlgorithm) = algo.params

function prepareSystemMatrix(reco::SinglePatchReconstructionParameter{L,S}) where {L<:AbstractSystemMatrixLoadingParameter, S<:AbstractLinearSolver}
  freqs, sf, grid = process(AbstractMPIRecoAlgorithm, reco.sfLoad, reco.sf, S, reco.arrayType)
  return freqs, sf, grid, reco.arrayType
end


AbstractImageReconstruction.take!(algo::SinglePatchReconstructionAlgorithm) = Base.take!(algo.output)

function process(algo::SinglePatchReconstructionAlgorithm, params::Union{A, ProcessResultCache{<:A}}, f::MPIFile, args...) where A <: AbstractMPIPreProcessingParameters
  result = process(typeof(algo), params, f, args...)
  if eltype(algo.S) != eltype(result)
    @warn "System matrix and measurement have different element data type. Mapping measurment data to system matrix element type."
    result = map(eltype(algo.S),result)
  end
  result = adapt(algo.arrayType, result)
  return result
end


function process(algo::SinglePatchReconstructionAlgorithm, params::SinglePatchReconstructionParameter, u)
  weights = adapt(algo.arrayType, process(algo, params.weightingParams, u))

  B = getLinearOperator(algo, params)

  solver = LeastSquaresParameters(op = B, S = algo.S, reg = params.reg, solverParams = params.solverParams, weights = weights)

  result = process(algo, solver, u)

  return gridresult(result, algo.grid, algo.sf)
end

function process(algo::SinglePatchReconstructionAlgorithm, params::Union{W, ProcessResultCache{W}}, u) where {W <: AbstractWeightingParameters} 
  result = process(typeof(algo), params, algo.freqs, algo.sf)
  if !isnothing(result)
    result = map(real(eltype(algo.S)), result)
  end
  return nothing
end

function getLinearOperator(algo::SinglePatchReconstructionAlgorithm, params::SinglePatchReconstructionParameter{<:DenseSystemMatixLoadingParameter, S}) where {S}
  return nothing
end

function getLinearOperator(algo::SinglePatchReconstructionAlgorithm, params::SinglePatchReconstructionParameter{<:SparseSystemMatrixLoadingParameter, S}) where {S}
  return process(algo, params.sfLoad, eltype(algo.S), algo.arrayType, tuple(shape(algo.grid)...))
end
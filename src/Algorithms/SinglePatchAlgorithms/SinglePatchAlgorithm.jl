Base.@kwdef struct SinglePatchReconstructionParameter{L<:AbstractSystemMatrixLoadingParameter, SL<:AbstractLinearSolver,
   arrT <: AbstractArray, SP<:AbstractSolverParameters{SL}, R<:AbstractRegularization, W<:AbstractWeightingParameters} <: AbstractSinglePatchReconstructionParameters
  # File
  sf::MPIFile
  sfLoad::Union{L, AbstractUtilityReconstructionParameters{L}}
  arrayType::Type{arrT} = Array
  # Solver
  solverParams::SP
  reg::Vector{R} = AbstractRegularization[]
  weightingParams::Union{W, AbstractUtilityReconstructionParameters{W}} = NoWeightingParameters()
end

@reconstruction constructor = false mutable struct SinglePatchReconstructionAlgorithm{P <: AbstractSinglePatchAlgorithmParameters, SM, arrT <: AbstractArray, vecT <: arrT} <: AbstractSinglePatchReconstructionAlgorithm
  @parameter params::P
  # Could also do reconstruction progress meter here
  sf::Union{MPIFile, Vector{MPIFile}}
  S::SM
  weights::Union{Nothing, vecT} = nothing
  arrayType::Type{arrT}
  grid::RegularGridPositions
  freqs::Vector{CartesianIndex{2}}
end

function SinglePatchReconstruction(params::SinglePatchParameters{<:AbstractMPIPreProcessingParameters, R, PT}) where {R<:AbstractSinglePatchReconstructionParameters, PT <:AbstractMPIPostProcessingParameters}
  return SinglePatchReconstructionAlgorithm(params)
end
function SinglePatchReconstructionAlgorithm(params::SinglePatchParameters{<:AbstractMPIPreProcessingParameters, R, PT}) where {R<:AbstractSinglePatchReconstructionParameters, PT <:AbstractMPIPostProcessingParameters}
  freqs, S, grid, arrayType = prepareSystemMatrix(params.reco)
  weights = prepareWeights(params.reco, freqs, S)
  return SinglePatchReconstructionAlgorithm{typeof(params), typeof(S), arrayType, typeof(arrayType{real(eltype(S))}(undef, 0))}(params, params.reco.sf, S, weights, arrayType, grid, freqs, @reconstruction_internals SinglePatchReconstructionAlgorithm)
end
recoAlgorithmTypes(::Type{SinglePatchReconstruction}) = SystemMatrixBasedAlgorithm()

function prepareSystemMatrix(reco::SinglePatchReconstructionParameter{L,S}) where {L<:AbstractSystemMatrixLoadingParameter, S<:AbstractLinearSolver}
  freqs, sf, grid = reco.sfLoad(AbstractMPIRecoAlgorithm, reco.sf, S, reco.arrayType)
  return freqs, sf, grid, reco.arrayType
end

function prepareWeights(reco::SinglePatchReconstructionParameter{L,S,arrT,SP,R,W}, freqs, sf) where {L, S, arrT, SP, R, W<:AbstractWeightingParameters}
  return reco.weightingParams(AbstractMPIRecoAlgorithm, freqs, sf, nothing, reco.arrayType)
end

function (params::SinglePatchParameters)(algo::SinglePatchReconstructionAlgorithm, data::MPIFile)
  return params(algo, data, algo.freqs)
end

function (params::Union{A, ProcessResultCache{<:A}})(algo::SinglePatchReconstructionAlgorithm, f::MPIFile, args...) where A <: AbstractMPIPreProcessingParameters
  result = params(typeof(algo), f, args...)
  if eltype(algo.S) != eltype(result)
    @warn "System matrix and measurement have different element data type. Mapping measurment data to system matrix element type."
    result = map(eltype(algo.S),result)
  end
  result = adapt(algo.arrayType, result)
  return result
end


function (params::SinglePatchReconstructionParameter)(algo::SinglePatchReconstructionAlgorithm, u)
  weights = params.weightingParams(algo, u, WeightingType(params.weightingParams))

  B = getLinearOperator(algo, params)

  solver = LeastSquaresParameters(op = B, S = algo.S, reg = params.reg, solverParams = params.solverParams, weights = weights)

  result = solver(algo, u)

  return gridresult(result, algo.grid, algo.sf)
end

function (params::Union{W, ProcessResultCache{W}})(algo::SinglePatchReconstructionAlgorithm, u, ::MeasurementBasedWeighting) where W<:AbstractWeightingParameters
  return params(typeof(algo), algo.freqs, algo.S, u, algo.arrayType)
end


function (params::Union{W, ProcessResultCache{W}})(algo::SinglePatchReconstructionAlgorithm, u, ::SystemMatrixBasedWeighting) where W<:AbstractWeightingParameters
  return algo.weights
end

function getLinearOperator(algo::SinglePatchReconstructionAlgorithm, params::SinglePatchReconstructionParameter)
  return nothing
end

function getLinearOperator(algo::SinglePatchReconstructionAlgorithm, params::SinglePatchReconstructionParameter{<:SparseSystemMatrixLoadingParameter})
  return params.sfLoad(algo, algo.sf, eltype(algo.S), algo.arrayType, tuple(shape(algo.grid)...))
end
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

Base.@kwdef mutable struct SinglePatchReconstructionAlgorithm{P, SM, arrT <: AbstractArray, vecT <: arrT} <: AbstractSinglePatchReconstructionAlgorithm where {P<:AbstractSinglePatchAlgorithmParameters}
  params::P
  # Could also do reconstruction progress meter here
  sf::Union{MPIFile, Vector{MPIFile}}
  S::SM
  weights::Union{Nothing, vecT} = nothing
  arrayType::Type{arrT}
  grid::RegularGridPositions
  freqs::Vector{CartesianIndex{2}}
  output::Channel{Any}
end

function SinglePatchReconstruction(params::SinglePatchParameters{<:AbstractMPIPreProcessingParameters, R, PT}) where {R<:AbstractSinglePatchReconstructionParameters, PT <:AbstractMPIPostProcessingParameters}
  return SinglePatchReconstructionAlgorithm(params)
end
function SinglePatchReconstructionAlgorithm(params::SinglePatchParameters{<:AbstractMPIPreProcessingParameters, R, PT}) where {R<:AbstractSinglePatchReconstructionParameters, PT <:AbstractMPIPostProcessingParameters}
  freqs, S, grid, arrayType = prepareSystemMatrix(params.reco)
  weights = prepareWeights(params.reco, freqs, S)
  return SinglePatchReconstructionAlgorithm{typeof(params), typeof(S), arrayType, typeof(arrayType{real(eltype(S))}(undef, 0))}(params, params.reco.sf, S, weights, arrayType, grid, freqs, Channel{Any}(Inf))
end
recoAlgorithmTypes(::Type{SinglePatchReconstruction}) = SystemMatrixBasedAlgorithm()
AbstractImageReconstruction.parameter(algo::SinglePatchReconstructionAlgorithm) = algo.params

function prepareSystemMatrix(reco::SinglePatchReconstructionParameter{L,S}) where {L<:AbstractSystemMatrixLoadingParameter, S<:AbstractLinearSolver}
  freqs, sf, grid = process(AbstractMPIRecoAlgorithm, reco.sfLoad, reco.sf, S, reco.arrayType)
  return freqs, sf, grid, reco.arrayType
end

function prepareWeights(reco::SinglePatchReconstructionParameter{L,S,arrT,SP,R,W}, freqs, S_) where {L, S, arrT, SP, R, W<:AbstractWeightingParameters}
  return process(AbstractMPIRecoAlgorithm, reco.weightingParams, freqs, reco.sf, S_, nothing, reco.arrayType)
end

Base.lock(algo::SinglePatchReconstructionAlgorithm) = lock(algo.output)
Base.unlock(algo::SinglePatchReconstructionAlgorithm) = unlock(algo.output)
Base.isready(algo::SinglePatchReconstructionAlgorithm) = isready(algo.output)
Base.wait(algo::SinglePatchReconstructionAlgorithm) = wait(algo.output)
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
  weights = process(algo, params.weightingParams, u, WeightingType(params.weightingParams))

  B = getLinearOperator(algo, params)

  solver = LeastSquaresParameters(op = B, S = algo.S, reg = params.reg, solverParams = params.solverParams, weights = weights)

  result = process(algo, solver, u)

  return gridresult(result, algo.grid, algo.sf)
end

function process(algo::SinglePatchReconstructionAlgorithm, params::Union{W, ProcessResultCache{W}}, u, ::MeasurementBasedWeighting) where W<:AbstractWeightingParameters
  return process(typeof(algo), params, algo.freqs, algo.S, u, algo.arrayType)
end


function process(algo::SinglePatchReconstructionAlgorithm, params::Union{W, ProcessResultCache{W}}, u, ::SystemMatrixBasedWeighting) where W<:AbstractWeightingParameters
  return algo.weights
end

function getLinearOperator(algo::SinglePatchReconstructionAlgorithm, params::SinglePatchReconstructionParameter{<:DenseSystemMatixLoadingParameter, S}) where {S}
  return nothing
end

function getLinearOperator(algo::SinglePatchReconstructionAlgorithm, params::SinglePatchReconstructionParameter{<:SparseSystemMatrixLoadingParameter, S}) where {S}
  return process(algo, params.sfLoad, eltype(algo.S), algo.arrayType, tuple(shape(algo.grid)...))
end
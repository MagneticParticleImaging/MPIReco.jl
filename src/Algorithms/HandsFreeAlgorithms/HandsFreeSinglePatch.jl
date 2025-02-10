export SinglePatchHandsFreeReconstructionParameter
Base.@kwdef struct SinglePatchHandsFreeReconstructionParameter{L<:AbstractSystemMatrixLoadingParameter,
  arrT <: AbstractArray, SP<:HandsFreeSolverParameters, W<:AbstractWeightingParameters} <: AbstractSinglePatchReconstructionParameters
 # File
 sf::MPIFile
 sfLoad::Union{L, ProcessResultCache{L}}
 arrayType::Type{arrT} = Array
 # Solver
 solverParams::SP
 #reg::Vector{R} = AbstractRegularization[]
 weightingParams::Union{W, ProcessResultCache{W}} = NoWeightingParameters()
end

function prepareSystemMatrix(reco::SinglePatchHandsFreeReconstructionParameter{L}) where {L<:AbstractSystemMatrixLoadingParameter}
  freqs, sf, grid = process(AbstractMPIRecoAlgorithm, reco.sfLoad, reco.sf, Kaczmarz, reco.arrayType)
  return freqs, sf, grid, reco.arrayType
end

function prepareWeights(reco::SinglePatchHandsFreeReconstructionParameter{L,arrT,SP,W}, freqs, sf) where {L, arrT, SP, W<:AbstractWeightingParameters}
  return process(AbstractMPIRecoAlgorithm, reco.weightingParams, freqs, sf, nothing, reco.arrayType)
end

function process(algo::SinglePatchReconstructionAlgorithm, params::SinglePatchHandsFreeReconstructionParameter, u)
  weights = process(algo, params.weightingParams, u, WeightingType(params.weightingParams))

  B = getLinearOperator(algo, params)

  solver = LeastSquaresParameters(op = B, S = algo.S, reg = L2Regularization[], solverParams = params.solverParams, weights = weights)

  snr = real(eltype(algo.S)).(vec(MPIFiles.getCalibSNR(algo.sf)[algo.freqs, 1]))

  result = process(algo, solver, u, snr)

  return gridresult(result, algo.grid, algo.sf)
end

function getLinearOperator(algo::SinglePatchReconstructionAlgorithm, params::SinglePatchHandsFreeReconstructionParameter{<:DenseSystemMatixLoadingParameter, S}) where {S}
  return nothing
end

function getLinearOperator(algo::SinglePatchReconstructionAlgorithm, params::SinglePatchHandsFreeReconstructionParameter{<:SparseSystemMatrixLoadingParameter, S}) where {S}
  return process(algo, params.sfLoad, eltype(algo.S), algo.arrayType, tuple(shape(algo.grid)...))
end
export SinglePatchHandsFreeReconstructionParameter
Base.@kwdef struct SinglePatchHandsFreeReconstructionParameter{L<:AbstractSystemMatrixLoadingParameter,
  arrT <: AbstractArray, SP<:HandsFreeSolverParameters, W<:AbstractWeightingParameters} <: AbstractSinglePatchReconstructionParameters
 # File
 sf::MPIFile
 sfLoad::Union{L, AbstractUtilityReconstructionParameters{L}}
 arrayType::Type{arrT} = Array
 # Solver
 solverParams::SP
 #reg::Vector{R} = AbstractRegularization[]
 weightingParams::Union{W, AbstractUtilityReconstructionParameters{W}} = NoWeightingParameters()
end

function prepareSystemMatrix(reco::SinglePatchHandsFreeReconstructionParameter{L}) where {L<:AbstractSystemMatrixLoadingParameter}
  freqs, sf, grid = reco.sfLoad(AbstractMPIRecoAlgorithm, reco.sf, Kaczmarz, reco.arrayType)
  return freqs, sf, grid, reco.arrayType
end

function prepareWeights(reco::SinglePatchHandsFreeReconstructionParameter{L,arrT,SP,W}, freqs, sf) where {L, arrT, SP, W<:AbstractWeightingParameters}
  return reco.weightingParams(AbstractMPIRecoAlgorithm, reco.weightingParams, freqs, sf, nothing, reco.arrayType)
end

function (params::SinglePatchHandsFreeReconstructionParameter)(algo::SinglePatchReconstructionAlgorithm, u)
  weights = params.weightingParams(algo, u, WeightingType(params.weightingParams))

  B = getLinearOperator(algo, params)

  solver = LeastSquaresParameters(op = B, S = algo.S, reg = L2Regularization[], solverParams = params.solverParams, weights = weights)

  snr = real(eltype(algo.S)).(vec(MPIFiles.getCalibSNR(algo.sf)[algo.freqs, 1]))

  result = solver(algo, u, snr)

  return gridresult(result, algo.grid, algo.sf)
end

function getLinearOperator(algo::SinglePatchReconstructionAlgorithm, params::SinglePatchHandsFreeReconstructionParameter)
  return nothing
end

function getLinearOperator(algo::SinglePatchReconstructionAlgorithm, params::SinglePatchHandsFreeReconstructionParameter{<:SparseSystemMatrixLoadingParameter})
  return params.sfLoad(algo, eltype(algo.S), algo.arrayType, tuple(shape(algo.grid)...))
end
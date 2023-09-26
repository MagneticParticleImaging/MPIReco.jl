export SinglePatchBGEstimationAlgorithm, SinglePatchBGEstimationReconstructionParameter
Base.@kwdef struct SinglePatchBGEstimationReconstructionParameter{L<:DenseSystemMatixLoadingParameter,
  SP<:AbstractSolverParameters} <: AbstractSinglePatchReconstructionParameters
  # File
  sf::MPIFile
  sfLoad::L
  # Solver
  solverParams::SP
  λ::Float32
  β::Float32
  # weightingType::WeightingType = WeightingType.None
  bgDict::BGDictParameter
end

Base.@kwdef mutable struct SinglePatchBGEstimationAlgorithm{P} <: AbstractSinglePatchReconstructionAlgorithm where {P<:AbstractSinglePatchAlgorithmParameters}
  params::P
  # Could also do reconstruction progress meter here
  sf::Union{MPIFile,Vector{MPIFile}}
  S::AbstractArray
  bgDict::AbstractArray
  grid::RegularGridPositions
  freqs::Vector{CartesianIndex{2}}
  output::Channel{Any}
end

function SinglePatchReconstruction(params::SinglePatchParameters{<:AbstractMPIPreProcessingParameters,<:SinglePatchBGEstimationReconstructionParameter,PT}) where {PT<:AbstractMPIPostProcessingParameters}
  return SinglePatchBGEstimationAlgorithm(params)
end
function SinglePatchBGEstimationAlgorithm(params::SinglePatchParameters{<:AbstractMPIPreProcessingParameters,R,PT}) where {R<:SinglePatchBGEstimationReconstructionParameter,PT<:AbstractMPIPostProcessingParameters}
  freqs, S, grid = prepareSystemMatrix(params.reco)
  return SinglePatchBGEstimationAlgorithm(params, params.reco.sf, S, process(SinglePatchBGEstimationAlgorithm, freqs, params.reco.bgDict), grid, freqs, Channel{Any}(Inf))
end
recoAlgorithmTypes(::Type{SinglePatchBGEstimationAlgorithm}) = SystemMatrixBasedAlgorithm()
AbstractImageReconstruction.parameter(algo::SinglePatchBGEstimationAlgorithm) = algo.origParam

function prepareSystemMatrix(reco::SinglePatchBGEstimationReconstructionParameter{L}) where {L<:AbstractSystemMatrixLoadingParameter}
  freqs, sf, grid = process(AbstractMPIRecoAlgorithm, reco.sfLoad, reco.sf)
  sf, grid = prepareSF(Kaczmarz, sf, grid)
  return freqs, sf, grid
end

AbstractImageReconstruction.take!(algo::SinglePatchBGEstimationAlgorithm) = Base.take!(algo.output)

function AbstractImageReconstruction.put!(algo::SinglePatchBGEstimationAlgorithm, data::MPIFile)
  consistenceCheck(algo.sf, data)

  result = process(algo, algo.params, data)

  # Create Image (maybe image parameter as post params?)
  # TODO make more generic to apply to other pre/reco params as well (pre.numAverage main issue atm)
  pixspacing = (spacing(algo.grid) ./ acqGradient(data)[1] .* acqGradient(algo.sf)[1]) * 1000u"mm"
  offset = (ffPos(data) .- 0.5 .* calibFov(algo.sf)) * 1000u"mm" .+ 0.5 .* pixspacing
  dt = acqNumAverages(data) * dfCycle(data) * algo.params.pre.numAverages * 1u"s"
  im = makeAxisArray(result, pixspacing, offset, dt)
  result = ImageMeta(im, generateHeaderDict(algo.sf, data))

  Base.put!(algo.output, result)
end


function process(algo::SinglePatchBGEstimationAlgorithm, params::AbstractMPIPreProcessingParameters, f::MPIFile)
  result = process(typeof(algo), params, f)
  if eltype(algo.S) != eltype(result)
    @warn "System matrix and measurement have different element data type. Mapping measurment data to system matrix element type."
    result = map(eltype(algo.S), result)
  end
  return result
end


function process(algo::SinglePatchBGEstimationAlgorithm, params::SinglePatchBGEstimationReconstructionParameter, u::Array)
  weights = nothing # getWeights(...)

  # Prepare Regularization
  reg = L2Regularization(Float32(params.λ))
  λ = RegularizedLeastSquares.normalize(params.solverParams.normalizeReg, reg, algo.S, nothing) # Scaling happens outside
  reg = L2Regularization(Float32(1))

  N = size(algo.S, 2)
  G = transpose(cat(Float32(1 / (sqrt(λ))) * transpose(algo.S), Float32(1 / (sqrt(params.β))) * transpose(algo.bgDict), dims=1))

  constraintMask = zeros(Bool, size(G, 2))
  constraintMask[1:N] .= 1
  # Enforce no normalization
  solverParams = fromKwargs(typeof(params.solverParams); toKwargs(params.solverParams, overwrite = Dict{Symbol, Any}(:normalizeReg => NoNormalization()))...)
  solverParams = ConstraintMaskedSolverParameters(;constraintMask = constraintMask, params = params.solverParams)

  solver = LeastSquaresParameters(solver = Kaczmarz, S = G, reg = [reg], solverParams = solverParams)

  temp = process(algo, solver, u)

  result = zeros(eltype(temp), N, size(temp, 2))

  for l = 1:size(temp, 2)
    d = temp[1:N, l] ./ sqrt(λ)
    result[:,l] = d
  end

  return gridresult(result, algo.grid, algo.sf)
end
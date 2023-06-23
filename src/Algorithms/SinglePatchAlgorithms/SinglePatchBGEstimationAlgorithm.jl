Base.@kwdef struct SinglePatchBGEstimationReconstructionParameter{L<:AbstractSystemMatrixLoadingParameter,S<:AbstractLinearSolver,
  SP<:AbstractSolverParameters,R<:AbstractRegularization} <: AbstractSinglePatchReconstructionParameters
  # File
  sf::MPIFile
  sfLoad::L
  # Solver
  solver::Type{S}
  solverParams::SP
  reg::Vector{R} = AbstractRegularization[]
  # weightingType::WeightingType = WeightingType.None
  bgDict::BGDictParameter
end

Base.@kwdef mutable struct SinglePatchBGEstimationAlgorithm{P} <: AbstractSinglePatchReconstructionAlgorithm where {P<:AbstractSinglePatchAlgorithmParameters}
  params::P
  # Could also do reconstruction progress meter here
  origParam::Union{AbstractSinglePatchAlgorithmParameters,Nothing} = nothing
  sf::Union{MPIFile,Vector{MPIFile}}
  S::AbstractArray
  bgDict::AbstractArray
  grid::RegularGridPositions
  freqs::Vector{Int64}
  output::Channel{Any}
end

function SinglePatchReconstruction(params::SinglePatchParameters{<:AbstractPreProcessingParameters,<:SinglePatchBGEstimationReconstructionParameter,PT}) where {PT<:AbstractPostProcessingParameters}
  return SinglePatchBGEstimationAlgorithm(params)
end
function SinglePatchBGEstimationAlgorithm(params::SinglePatchParameters{<:AbstractPreProcessingParameters,R,PT}) where {R<:SinglePatchBGEstimationReconstructionParameter,PT<:AbstractPostProcessingParameters}
  freqs, S, grid = prepareSystemMatrix(params.pre, params.reco)
  filter = fromKwargs(FrequencyFilteredPreProcessingParameters; frequencies=freqs, toKwargs(params.pre; flatten=DataType[])...)
  filteredParams = SinglePatchParameters(filter, params.reco, params.post)
  return SinglePatchBGEstimationAlgorithm(filteredParams, params, params.reco.sf, S, process(SinglePatchBGEstimationAlgorithm, freqs, params.reco.bgDict), grid, freqs, Channel{Any}(Inf))
end
recoAlgorithmTypes(::Type{SinglePatchReconstruction}) = SystemMatrixBasedAlgorithm()
RecoUtils.parameter(algo::SinglePatchBGEstimationAlgorithm) = algo.origParam

function prepareSystemMatrix(pre::AbstractPreProcessingParameters, reco::SinglePatchBGEstimationReconstructionParameter{L,S}) where {L<:AbstractSystemMatrixLoadingParameter,S<:AbstractLinearSolver}
  params = fromKwargs(PreProcessedSystemMatrixLoadingParameter; pre=pre, sm=reco.sfLoad)
  freqs, sf, grid = process(AbstractMPIReconstructionAlgorithm, reco.sf, params)
  sf, grid = prepareSF(S, sf, grid)
  return freqs, sf, grid
end

RecoUtils.take!(algo::SinglePatchBGEstimationAlgorithm) = Base.take!(algo.output)

function RecoUtils.put!(algo::SinglePatchBGEstimationAlgorithm, data::MPIFile)
  consistenceCheck(algo.sf, data)

  result = process(algo, data, algo.params)

  # Create Image (maybe image parameter as post params?)
  # TODO make more generic to apply to other pre/reco params as well (pre.numAverage main issue atm)
  pixspacing = (spacing(algo.grid) ./ acqGradient(data)[1] .* acqGradient(algo.sf)[1]) * 1000u"mm"
  offset = (ffPos(data) .- 0.5 .* calibFov(algo.sf)) * 1000u"mm" .+ 0.5 .* pixspacing
  dt = acqNumAverages(data) * dfCycle(data) * algo.params.pre.numAverages * 1u"s"
  im = makeAxisArray(result, pixspacing, offset, dt)
  result = ImageMeta(im, generateHeaderDict(algo.sf, data))

  Base.put!(algo.output, result)
end

function RecoUtils.similar(algo::SinglePatchBGEstimationAlgorithm, data::MPIFile)
  # TODO Check num receive channel

  # This isnt a nice structure atm, because the corrections of the params cant be done individually
  # and in certain parts should only be applied after construction

  # Check if freq still valid
  freqs = algo.freqs
  S = algo.S
  grid = algo.grid
  params = parameter(algo)
  if rxNumSamplingPoints(algo.sf) == rxNumSamplingPoints(data)
    # Ensure that no frequencies are used that are not present in the measurement
    freqParams = fromKwargs(PreProcessedFrequencyFilterParameter; toKwargs([params.pre, params.reco]; flatten=DataType[AbstractSystemMatrixLoadingParameter])...)
    measFreqs = process(AbstractMPIReconstructionAlgorithm, data, freqParams)
    freqs = intersect(algo.freqs, measFreqs)
    if freqs != algo.freqs
      S, grid = getSF(params.reco.sf, freqs, nothing; toKwargs(params.pre)...)
      S, grid = prepareSF(params.reco.solver, S, grid)
    end
  end

  numPeriodGrouping = params.pre.numPeriodGrouping
  if rxNumSamplingPoints(params.reco.sf) > rxNumSamplingPoints(data)
    numPeriodGrouping = rxNumSamplingPoints(params.reco.sf) ÷ rxNumSamplingPoints(data)
  end
  numPeriodAverages = params.pre.numPeriodAverages
  if acqNumPeriodsPerFrame(params.reco.sf) < acqNumPeriodsPerFrame(data)
    numPeriodAverages = acqNumPeriodsPerFrame(data) ÷ (acqNumPeriodsPerFrame(params.reco.sf) * numPeriodGrouping)
  end

  pre = fromKwargs(FrequencyFilteredPreProcessingParameters; toKwargs(params.pre, flatten=DataType[], overwrite=Dict{Symbol,Any}(:numPeriodGrouping => numPeriodGrouping, :numPeriodAverages => numPeriodAverages, :frequencies => freqs))...)
  reco = RecoUtils.similar(algo, data, algo.params.reco)
  post = RecoUtils.similar(algo, data, algo.params.post)

  params = SinglePatchParameters(pre, reco, post)

  result = SinglePatchBGEstimationAlgorithm(params, params, algo.params.reco.sf, S, grid, freqs, Channel{Any}(Inf))

  return result
end

function RecoUtils.process(algo::SinglePatchBGEstimationAlgorithm, f::MPIFile, params::AbstractPreProcessingParameters)
  result = process(typeof(algo), f, params)
  if eltype(algo.S) != eltype(result)
    @warn "System matrix and measurement have different element data type. Mapping measurment data to system matrix element type."
    result = map(eltype(algo.S), result)
  end
  return result
end


function RecoUtils.process(algo::SinglePatchBGEstimationAlgorithm, u::Array, params::SinglePatchReconstructionParameter)
  weights = nothing # getWeights(...)

  B = getLinearOperator(algo, params)

  G = transpose(cat(Float32(1 / (sqrt(λ))) * transpose(S), Float32(1 / (sqrt(β))) * transpose(bgDict), dims=1))

  solver = LeastSquaresParameters(params.solver, nothing, G, params.reg, params.solverParams)

  temp = process(algo, u, solver)

  N = size(algo.S, 2)
  result = zeros(eltype(temp), N, size(temp, 2))

  for l = 1:size(temp, 2)
    d = temp[1:N, l] ./ sqrt(λ)
    if !isnothing(B)
      d[:] = B*d
    end
    result[:,l] = d
  end

  numcolors = 1
  if isa(algo.sf, AbstractVector) || isa(algo.sf, MultiContrastFile)
    numcolors = length(algo.sf)
  end
  shp = shape(algo.grid)
  cArray = Array{Float32}(undef, numcolors, shp..., size(result)[end])
  result = reshape(result, reduce(*, shp), numcolors, :)
  result = permutedims(result, [2, 1, 3])
  cArray[:] = result[:]

  return cArray
end

function getLinearOperator(algo::SinglePatchBGEstimationAlgorithm, params::SinglePatchReconstructionParameter{<:DenseSystemMatixLoadingParameter,S}) where {S}
  return linearOperator(nothing, shape(algo.grid), eltype(algo.S))
end

function getLinearOperator(algo::SinglePatchBGEstimationAlgorithm, params::SinglePatchReconstructionParameter{<:SparseSystemMatrixLoadingParameter,S}) where {S}
  return linearOperator(params.sfLoad.sparseTrafo, shape(algo.grid), eltype(algo.S))
end
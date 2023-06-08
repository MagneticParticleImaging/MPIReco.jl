Base.@kwdef struct SinglePatchReconstructionParameter{L<:AbstractSystemMatrixLoadingParameter, S<:AbstractLinearSolver, SP<:AbstractSolverIterationParameters, R<:AbstractRegularization} <: AbstractSinglePatchReconstructionParameters
  # File
  sf::MPIFile
  sfLoad::L
  # Solver
  solver::Type{S}
  solverParams::SP
  reg::Vector{R} = AbstractRegularization[]
  # weightingType::WeightingType = WeightingType.None
end

Base.@kwdef mutable struct SinglePatchParameters{PR<:AbstractPreProcessingParameters, R<:AbstractSinglePatchReconstructionParameters, PT<:AbstractPostProcessingParameters} <: AbstractRecoAlgorithmParameters
  pre::Union{PR}
  reco::Union{R}
  post::Union{PT} = NoPostProcessing() 
end

Base.@kwdef mutable struct SinglePatchReconstructionAlgorithm{PR, R, PT} <: AbstractMPIReconstructionAlgorithm where {PR, R, PT<:AbstractMPIRecoParameters}
  params::SinglePatchParameters{PR, R, PT}
  # Could also do reconstruction progress meter here
  sf::Union{MPIFile, Vector{MPIFile}}
  S::AbstractArray
  grid::RegularGridPositions
  freqs::Vector{Int64}
  output::Channel{Any}
end

function SinglePatchReconstruction(params::SinglePatchParameters{<:CommonPreProcessingParameters, R, PT}) where {R<:AbstractSinglePatchReconstructionParameters, PT <:AbstractPostProcessingParameters}
  return SinglePatchReconstructionAlgorithm(params)
end
function SinglePatchReconstructionAlgorithm(params::SinglePatchParameters{<:CommonPreProcessingParameters, R, PT}) where {R<:AbstractSinglePatchReconstructionParameters, PT <:AbstractPostProcessingParameters}
  freqs, S, grid = prepareSystemMatrix(params.pre, params.reco)
  filter = fromKwargs(FrequencyFilteredPreProcessingParameters; frequencies = freqs, toKwargs(params.pre; flatten = DataType[])...)
  filteredParams = SinglePatchParameters(filter, params.reco, params.post)
  return SinglePatchReconstructionAlgorithm(filteredParams, params.reco.sf, S, grid, freqs, Channel{Any}(Inf))
end
recoAlgorithmTypes(::Type{SinglePatchReconstruction}) = SystemMatrixBasedAlgorithm()

function prepareSystemMatrix(pre::CommonPreProcessingParameters, reco::SinglePatchReconstructionParameter{L,S}) where {L<:AbstractSystemMatrixLoadingParameter, S<:AbstractLinearSolver}
  params = fromKwargs(PreProcessedSystemMatrixLoadingParameter; toKwargs([pre, reco])..., sm = reco.sfLoad)
  freqs, sf, grid = process(AbstractMPIReconstructionAlgorithm, reco.sf, params)
  sf, grid = prepareSF(S, sf, grid) 
  return freqs, sf, grid
end


RecoUtils.take!(algo::SinglePatchReconstructionAlgorithm) = Base.take!(algo.output)

function RecoUtils.put!(algo::SinglePatchReconstructionAlgorithm, data::MPIFile)
  consistenceCheck(algo.params.reco.sf, data)
  
  result = process(algo, data, algo.params.pre)
  
  result = process(algo, result, algo.params.reco)
  
  result = process(algo, result, algo.params.post)
  
  # Create Image (maybe image parameter as post params?)
  # TODO make more generic to apply to other pre/reco params as well (pre.numAverage main issue atm)
  pixspacing = (spacing(algo.grid) ./ acqGradient(data)[1] .* acqGradient(algo.sf)[1])*1000u"mm"
  offset = (ffPos(data) .- 0.5 .* calibFov(algo.sf))*1000u"mm" .+ 0.5 .* pixspacing
  dt = acqNumAverages(data)*dfCycle(data)*algo.params.pre.numAverages*1u"s"
  im = makeAxisArray(result, pixspacing, offset, dt)
  result = ImageMeta(im, generateHeaderDict(algo.sf, data))

  Base.put!(algo.output, result)
end

function RecoUtils.similar(algo::SinglePatchReconstructionAlgorithm, data::MPIFile)
  # TODO Check num receive channel

  # This isnt a nice structure atm, because the corrections of the params cant be done individually
  # and in certain parts should only be applied after construction

  # Check if freq still valid
  freqs = algo.freqs
  S = algo.S
  grid = algo.grid
  if rxNumSamplingPoints(algo.params.reco.sf) == rxNumSamplingPoints(data)
    # Ensure that no frequencies are used that are not present in the measurement
    freqParams = fromKwargs(PreProcessedFrequencyFilterParameter; toKwargs([algo.params.pre, algo.params.reco]; flatten = DataType[AbstractSystemMatrixLoadingParameter])...)
    measFreqs = process(AbstractMPIReconstructionAlgorithm, data, freqParams)
    freqs = intersect(algo.freqs, measFreqs)
    if freqs != algo.freqs
      S, grid = getSF(algo.params.reco.sf, freqs, nothing; toKwargs(algo.params.pre)...)
      S, grid = prepareSF(algo.params.reco.solver, S, grid)
    end
  end

  numPeriodGrouping = algo.params.pre.numPeriodGrouping
  if rxNumSamplingPoints(algo.params.reco.sf) > rxNumSamplingPoints(data)
    numPeriodGrouping = rxNumSamplingPoints(algo.params.reco.sf) ÷ rxNumSamplingPoints(data)
  end
  numPeriodAverages = algo.params.pre.numPeriodAverages
  if acqNumPeriodsPerFrame(algo.params.reco.sf) < acqNumPeriodsPerFrame(data)
    numPeriodAverages = acqNumPeriodsPerFrame(data) ÷ (acqNumPeriodsPerFrame(algo.params.reco.sf) * numPeriodGrouping)
  end

  pre = fromKwargs(FrequencyFilteredPreProcessingParameters; toKwargs(algo.params.pre, overwrite = Dict{Symbol, Any}(:numPeriodGrouping => numPeriodGrouping, :numPeriodAverages => numPeriodAverages, :frequencies => freqs))...)
  reco = RecoUtils.similar(algo, data, algo.params.reco)
  post = RecoUtils.similar(algo, data, algo.params.post)

  result = SinglePatchReconstructionAlgorithm(SinglePatchParameters(pre, reco, post), S, grid, freqs, Channel{Any}(32))

  return result
end

function RecoUtils.process(algo::SinglePatchReconstructionAlgorithm, f::MPIFile, params::AbstractPreProcessingParameters)
  result = process(AbstractMPIReconstructionAlgorithm, f, params)
  if eltype(algo.S) != eltype(result)
    @warn "System matrix and measurement have different element data type. Mapping measurment data to system matrix element type."
    result = map(eltype(algo.S),result)
  end
  return result
end


function RecoUtils.process(algo::SinglePatchReconstructionAlgorithm, u::Array, params::SinglePatchReconstructionParameter)
  weights = nothing # getWeights(...)

  B = getLinearOperator(algo, params)

  solver = LeastSquaresParameters(params.solver, B, algo.S, params.reg, params.solverParams)

  result = process(AbstractMPIReconstructionAlgorithm, u, solver)

  numcolors = 1
  if isa(algo.sf,AbstractVector) || isa(algo.sf,MultiContrastFile)
    numcolors = length(algo.sf)
  end
  shp = shape(algo.grid)
  cArray = Array{Float32}(undef, numcolors, shp..., size(result)[end])
  result = reshape(result, reduce(*, shp),numcolors,:)
  result = permutedims(result, [2,1,3])
  cArray[:] = result[:]

  return cArray
end

function getLinearOperator(algo::SinglePatchReconstructionAlgorithm, params::SinglePatchReconstructionParameter{<:DenseSystemMatixLoadingParameter, S}) where {S}
  return linearOperator(nothing, shape(algo.grid), eltype(algo.S))
end

function getLinearOperator(algo::SinglePatchReconstructionAlgorithm, params::SinglePatchReconstructionParameter{<:SparseSystemMatrixLoadingParameter, S}) where {S}
  return linearOperator(params.sfLoad.sparseTrafo, shape(algo.grid), eltype(algo.S))
end
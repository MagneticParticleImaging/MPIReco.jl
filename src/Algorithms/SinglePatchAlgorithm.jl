Base.@kwdef struct SinglePatchReconstructionParameter{S<:AbstractLinearSolver} <: AbstractReconstructionParameters
  # File
  sf::MPIFile
  sparseTrafo::Union{Nothing, String} = nothing
  #saveTrafo::Bool = false
  recChannels::UnitRange{Int64} = 1:rxNumChannels(sf)
  # Freqs
  minFreq::Float64 = 0.0
  maxFreq::Float64 = rxBandwidth(sf)
  numUsedFreqs::Int64 = -1
  # SNR
  threshSNR::Float64=-1.0
  sortBySNR::Bool = false
  # Solver
  solver::Type{S}
  iterations::Int64=10
  enforceReal=false
  enforcePositive=true
  λ::Vector{Float64}=[0.1, 0.0]
  relativeLambda::Bool=true
  # Grid
  gridsize::Vector{Int64} = gridSizeCommon(sf)
  fov::Vector{Float64} = calibFov(sf)
  center::Vector{Float64} = [0.0,0.0,0.0]
  useDFFoV::Bool = false
  # Misc.
  deadPixels::Vector{Int64} = Int64[]
end

Base.@kwdef mutable struct SinglePatchParameters{PR<:AbstractPreProcessingParameters, R<:AbstractReconstructionParameters, PT<:AbstractPostProcessingParameters} <: AbstractRecoAlgorithmParameters
  pre::Union{PR,Nothing} = nothing
  reco::Union{R,Nothing} = nothing
  post::Union{PT,Nothing} = nothing
end

Base.@kwdef mutable struct SinglePatchReconstruction{PR, R, PT} <: AbstractMPIReconstructionAlgorithm where {PR, R, PT<:MPIRecoParameters}
  params::SinglePatchParameters{PR, R, PT}
  # Could also do reconstruction progress meter here
  S
  grid
  freqs
  output::Channel{Any}
end

function SinglePatchReconstruction(params::SinglePatchParameters)
  # Prepare system matrix based on pre and reco params
  freqs = getFreqs(params.pre, params.reco)
  S, grid = getSF(params.pre, params.reco, freqs)
  filter = FrequencyFilteredPreProcessingParameters(;frequencies = freqs, toKwargs(params.pre)...)
  filteredParams = SinglePatchParameters(filter, params.reco, params.post)
  return SinglePatchReconstruction(filteredParams, S, grid, freqs, Channel{Any}(32))
end

recoAlgorithmTypes(::Type{SinglePatchReconstruction}) = SystemMatrixBasedAlgorithm()

getFreqs(pre::AbstractPreProcessingParameters, reco::SinglePatchReconstructionParameter) = getFreqs(reco.sf, pre, reco)
function getFreqs(data::MPIFile, pre::AbstractPreProcessingParameters, reco::SinglePatchReconstructionParameter) 
  dict = toKwargs([pre, reco])
  # Hacky version, make nicer later (maybe change in MPIFiles)
  # TODO fix this
  return filterFrequencies(data; SNRThresh = dict[:threshSNR], minFreq = dict[:minFreq], maxFreq = dict[:maxFreq],
   recChannels = dict[:recChannels], sortBySNR = dict[:sortBySNR], numUsedFreqs = dict[:numUsedFreqs], 
   numPeriodAverages = dict[:numPeriodAverages], numPeriodGrouping = dict[:numPeriodGrouping])
end

# TODO Maybe do this on Type{T}
function getSF(pre::AbstractPreProcessingParameters, reco::SinglePatchReconstructionParameter, freqs)
  # TODO make this work with correct bg correction
  return getSF(reco.sf, freqs, reco.sparseTrafo, reco.solver; toKwargs([pre, reco]; ignore = [:bgCorrection])...)
end

RecoUtils.take!(algo::SinglePatchReconstruction) = Base.take!(algo.output)

function RecoUtils.put!(algo::SinglePatchReconstruction, data::MPIFile)
  consistenceCheck(algo.params.reco.sf, data)
  
  pre = process(algo, data, algo.params.pre)
  
  result = process(algo, pre, algo.params.reco)
  
  result = process(algo, result, algo.params.post)
  
  # Create Image (maybe image parameter as post params?)
  # TODO make more generic to apply to other pre/reco params as well (pre.numAverage main issue atm)
  pixspacing = (spacing(algo.grid) ./ acqGradient(data)[1] .* acqGradient(algo.params.reco.sf)[1])*1000u"mm"
  offset = (ffPos(data) .- 0.5 .* calibFov(algo.params.reco.sf))*1000u"mm" .+ 0.5 .* pixspacing
  dt = acqNumAverages(data)*dfCycle(data)*algo.params.pre.numAverages*1u"s"
  axis = ImageAxisParameter(pixspacing=pixspacing, offset=offset, dt = dt)
  imMeta = ImageMetadataSystemMatrixParameter(data, algo.params.reco.sf, algo.grid, axis)
  result = process(AbstractMPIReconstructionAlgorithm, result, imMeta)

  Base.put!(algo.output, result)
end

function RecoUtils.similar(algo::SinglePatchReconstruction, data::MPIFile)
  # TODO Check num receive channel

  # This isnt a nice structure atm, because the corrections of the params cant be done individually
  # and in certain parts should only be applied after construction

  # Check if freq still valid
  freqs = algo.freqs
  S = algo.S
  grid = algo.grid
  if rxNumSamplingPoints(algo.params.reco.sf) == rxNumSamplingPoints(data)
    # Ensure that no frequencies are used that are not present in the measurement
    freqs = intersect(algo.freqs, getFreqs(data, algo.params.pre, algo.params.reco))
    if freqs != algo.freqs
      S, grid = getSF(params.pre, params.reco, freqs)
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

  result = SinglePatchReconstruction(SinglePatchParameters(pre, reco, post), S, grid, freqs, Channel{Any}(32))

  return result
end

function RecoUtils.process(algo::SinglePatchReconstruction, f::MPIFile, params::AbstractPreProcessingParameters)
  result = process(AbstractMPIReconstructionAlgorithm, f, params)
  if eltype(algo.S) != eltype(result)
    @warn "System matrix and measurement have different element data type. Mapping measurment data to system matrix element type."
    result = map(eltype(algo.S),result)
  end
  return result
end


function RecoUtils.process(algo::SinglePatchReconstruction, u::Array, params::SinglePatchReconstructionParameter)
  weights = nothing # getWeights(...)

  B = linearOperator(params.sparseTrafo, shape(algo.grid), eltype(algo.S))

  N = size(algo.S, 2)
  M = div(length(algo.S), N)
  L = size(u)[end]
  u = reshape(u, M, L)
  c = zeros(N, L)

  λ = params.λ

  if sum(abs.(λ)) > 0 && params.solver != FusedLasso && params.relativeLambda
    trace = calculateTraceOfNormalMatrix(algo.S,weights)
    if isa(λ,AbstractVector) 
      λ[1:1] *= trace / N
    else
      λ *= trace / N
    end
    #setlambda(S,λ) dead code?
  end

  args = toKwargs(params)
  args[:λ] = λ
  args[:sparseTrafo] = B
  solv = createLinearSolver(params.solver, algo.S; weights=weights, args...)

  for l=1:L
    d = solve(solv, u[:, l])
    if !isnothing(B)
      d[:] = B*d
    end
    c[:, l] = real(d)
  end

  return c
end
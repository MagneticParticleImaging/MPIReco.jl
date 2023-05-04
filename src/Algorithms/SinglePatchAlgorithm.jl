Base.@kwdef struct SinglePatchReconstructionParameter{S<:AbstractLinearSolver} <: AbstractReconstructionParameters
  # File
  sf::MPIFile
  sparseTrafo::Union{Nothing, String} = nothing
  #saveTrafo::Bool = false
  recChannels::UnitRange{Int64} = 1:rxNumChannels(sf)
  # Freqs
  minFreq::Int64 = 0
  maxFreq::Int64 = rxBandwidth(sf)
  numUsedFreqs::Int64 = -1
  # SNR
  threshSNR::Float64=-1
  sortBySNR::Bool = false
  # Solver
  solver::Type{S}
  iterations::Int64=10
  enforceReal=true
  enforcePositive=true
  λ::Vector{Float64}=[0.1, 0.0]
  relativeLambda::Bool=true
  # Grid
  gridsize::Vector{Int64} = gridSizeCommon(bSF)
  fov::Vector{Float64} = calibFov(bSF)
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

Base.@kwdef mutable struct SinglePatchReconstruction <: AbstractMPIReconstructionAlgorithm
  params::SinglePatchParameters
  # Could also do reconstruction progress meter here
  S
  grid
  freqs
  output::Channel{Any}
end

function SinglePatchReconstruction(params::SinglePatchReconstructionParameter)
  # Prepare system matrix based on pre and reco params
  freqs = getFreqs(params.pre, params.reco)
  S, grid = getSF(params.pre, params.reco, freqs)
  # params.pre = FrequencyFilteredPreProcessingParameters(;freqs = freqs, toKwargs(params.pre))# TODO: define FrequencyFilteredPreProcessingParameters
  # At the end we want to be in a state where we can handle a series of put! requests
  return SinglePatchReconstruction(params, S, grid, freqs, Channel{Any}(32))
end

recoAlgorithmTypes(::Type{SinglePatchReconstruction}) = SystemMatrixBasedAlgorithm()

getFreqs(pre::CommonPreProcessingParameters, reco::SinglePatchReconstructionParameter) = filterFrequences(reco.sf, toKargs([pre, reco])...)

# TODO Maybe do this on Type{T}
function getSF(pre::CommonPreProcessingParameters, reco::SinglePatchReconstructionParameter, freqs)
  return getSF(reco.sf, freqs, reco.sparseTrafo, reco.solver; toKwargs([pre, reco])...)
end

take!(algo::SinglePatchReconstruction) = take!(algo.output)

function put!(algo::SinglePatchReconstruction, data)
  consistenceCheck(bSF, bMeas)
  pre = process(algo, data, algo.params.pre)
  result = process(algo, pre, algo.params.reco)
  if !isnothing(algo.params.post)
    result = process(algo, result, algo.params.post)
  end
  put!(algo.output, result)
end

function similar(algo::SinglePatchReconstruction, data::MPIFile)
  pre = similar(algo, data, algo.param.pre)
  reco = similar(algo, data, algo.param.reco)
  post = similar(algo, data, algo.parms.post)

  # TODO Check num receive channel

  # This isnt a nice structure atm, because the corrections of the params cant be done individually
  # and in certain parts should only be applied after construction

  # Check if freq still valid
  freqs = algo.freqs
  S = algo.S
  grid = algo.grid
  if rxNumSamplingPoints(reco.sf) == rxNumSamplingPoints(data)
    # Ensure that no frequencies are used that are not present in the measurement
    freqs = intersect(algo.freqs, filterFrequencies(data, toKwargs(pre)))
    if freqs != algo.freqs
      S, grid = getSF(params.pre, params.reco, freqs)
    end
  end

  result = SinglePatchReconstruction(SinglePatchReconstructionParameter(pre, reco, post), S, grid, freqs, Channel{Any}(32))

  # TODO this only works for common preprocessing, because i directly access fields
  # If S is processed and fits not to the measurements because of numPeriodsGrouping
  # or numPeriodAverages being applied we need to set these so that the 
  # measurements are loaded correctly
  if rxNumSamplingPoints(bSF) > rxNumSamplingPoints(bMeas)
    pre.numPeriodGrouping = rxNumSamplingPoints(bSF) ÷ rxNumSamplingPoints(bMeas)
  end
  if acqNumPeriodsPerFrame(bSF) < acqNumPeriodsPerFrame(bMeas)
    pre.numPeriodAverages = acqNumPeriodsPerFrame(bMeas) ÷ (acqNumPeriodsPerFrame(bSF) * numPeriodGrouping)
  end
  return result
end

function process(algo::SinglePatchReconstruction, f::MPIFile, params::AbstractPreProcessingParameters)
  result = process(AbstractMPIReconstructionAlgorithm, f, params)
  if eltype(algo.S) != eltype(result)
    @warn "System matrix and measurement have different element data type. Mapping measurment data to system matrix element type."
    result = map(eltype(algo.S),result)
  end
  return result
end


function process(algo::SinglePatchReconstruction, u::Array, params::SinglePatchReconstructionParameter)
  weights = nothing # getWeights(...)

  B = linearOperator(params.sparseTrafo, shape(grid), eltype(S))

  N = size(algo.S, 2)
  M = div(length(S), N)
  L = size(u)[end]
  u = reshape(u, M, L)
  c = zeros(N, L)

  λ = params.λ

  if sum(abs.(λ)) > 0 && solver != FusedLasso && relativeLambda
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
  solv = createLinearSolver(params.solver, S; weights=weights, args...)

  for l=1:L
    d = solve(solv, u[:, l])
    if !isnothing(B)
      d[:] = B*d
    end
    c[:, l] = real(d)
  end

  return c
end
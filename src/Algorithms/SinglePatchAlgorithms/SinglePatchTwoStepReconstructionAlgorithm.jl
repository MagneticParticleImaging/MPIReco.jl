Base.@kwdef struct SinglePatchTwoStepReconstructionParameters{S<:AbstractLinearSolver} <: AbstractSinglePatchReconstructionParameters
  # Threshhold
  Γ::Float64
  # File
  sf::MPIFile
  sparseTrafo::Union{Nothing, String} = nothing # TODO concrete options here?
  #saveTrafo::Bool = false
  recChannels::UnitRange{Int64} = 1:rxNumChannels(sf)
  # Freqs
  minFreq::Float64 = 0.0
  maxFreq::Float64 = rxBandwidth(sf)
  # SNR
  threshSNR_high::Float64=-1.0
  threshSNR_low::Float64=-1.0
  sortBySNR::Bool = false
  # Solver
  solver::Type{S}
  # TODO Maybe nest these too
  iterations_high::Int64=10
  iterations_low::Int64=3
  enforceReal=false
  enforcePositive=true
  λ_high::Vector{Float64}=[0.1]
  λ_low::Vector{Float64}=[0.1]
  relativeLambda::Bool=true
  # weightingType::WeightingType = WeightingType.None
  # Grid
  gridsize::Vector{Int64} = gridSizeCommon(sf)
  fov::Vector{Float64} = calibFov(sf)
  center::Vector{Float64} = [0.0,0.0,0.0]
  useDFFoV::Bool = false
  # Misc.
  deadPixels::Vector{Int64} = Int64[]  
end
Base.@kwdef mutable struct SinglePatchTwoStepReconstructionAlgorithm{PR, R, PT} <: AbstractMPIReconstructionAlgorithm where {PR, R, PT<:AbstractMPIRecoParameters}
  params::SinglePatchParameters{PR, R, PT}
  # Could also do reconstruction progress meter here
  algoHigh::SinglePatchReconstructionAlgorithm
  algoLow::SinglePatchReconstructionAlgorithm
  output::Channel{Any}
end


function SinglePatchReconstruction(params::SinglePatchParameters{PR, SinglePatchTwoStepReconstructionParameters{solv}, PT}) where {solv, PR<:AbstractPreProcessingParameters, PT <:AbstractPostProcessingParameters}
  recoHigh = fromKwargs(SinglePatchReconstructionParameter, overwrite = Dict{Symbol, Any}(:iterations => params.reco.iterations_high, :threshSNR => params.reco.threshSNR_high))
  recoLow = fromKwargs(SinglePatchReconstructionParameter, overwrite = Dict{Symbol, Any}(:iterations => params.reco.iterations_low, :threshSNR => params.reco.threshSNR_low))
  algoHigh = SinglePatchReconstruction(SinglePatchReconstructionParameter(params.pre, recoHigh, params.post))
  algoLow = SinglePatchReconstruction(SinglePatchReconstructionParameter(params.pre, recoLow, params.post))
  return SinglePatchTwoStepReconstructionAlgorithm(params, algoHigh, algoLow, Channel{Any}(Inf))
end

function RecoUtils.put!(algo::SinglePatchTwoStepReconstructionParameters, data::MPIFile)
  # First reco
  cPre = reconstruct(algo.algoHigh, data)

  # Thresholding
  cThresh = copy(cPre)
  cThresh[ abs.(cPre).< maximum(abs.(cPre))*algo.params.reco.Γ ] .= 0

  # Projection into raw data space
  uProj = map(ComplexF32,algoLow.S*vec(cThresh))

  # Subtraction
  uMeas_low = process(algo.algoLow, data, algo.algoLow.params.pre)
  uCorr = uMeas_low - uProj

  # Second reconstruction
  cPost = reconstruct(algo.algoLow, uCorr)

  # Addition
  result = cPost + cThresh

  Base.put!(algo.output, result)
end

function RecoUtils.similar(algo::SinglePatchTwoStepReconstructionAlgorithm, data::MPIFile)
  algoHigh = similar(algo.algoHigh, data)
  pre = fromKwargs(CommonPreProcessingParameters, algoHigh.params.pre)
  reco = RecoUtils.similar(algo, data, algo.params.reco)
  post = RecoUtils.similar(algo, data, algo.params.post)
  return SinglePatchReconstruction(SinglePatchReconstructionParameters(pre, reco, post))
end


### Alternative
#=

# Extend structs to accommodate for two step, multi color, single and more by allowing multiple freqs and SMs 
Base.@kwdef mutable struct SinglePatchReconstructionAlgorithm{PR, R, PT} <: AbstractMPIReconstructionAlgorithm where {PR, R, PT<:AbstractMPIRecoParameters}
  params::SinglePatchParameters{PR, R, PT}
  # Could also do reconstruction progress meter here
  S::Vector{AbstractArray}
  grid::RegularGridPositions
  freqs::Vector{Vector{Int64}}
  output::Channel{Any}
end

# In constructor and in put! dispatch on parametric type of reco args
function Base.put!(algo::SinglePatchAlgorithm{PR, SinglePatchTwoStepReconstructionParameters{S}, PT}) where {PR, S, PT}
  consistenceCheck(algo.params.reco.sf, data)
  pre = process(algo, data, algo.params.pre)

  # Step 1
  result = process(algo, pre, ...) # Build other reco params here
  result = process(algo, result, algo.params.post)
  ...  
  imMeta = ImageMetadataSystemMatrixParameter(data, algo.params.reco.sf, algo.grid, axis)
  result = process(AbstractMPIReconstructionAlgorithm, result, imMeta)

  # Step 4
  result = process(algo, pre, ...) # Build other reco params for second reco here
  result = process(algo, result, algo.params.post)
  ...

    # Addition
  result = cPost + cThresh

  Base.put!(algo.output, result)
end


=#
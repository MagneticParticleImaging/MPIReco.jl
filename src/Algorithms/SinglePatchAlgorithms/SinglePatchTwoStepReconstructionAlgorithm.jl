export SinglePatchTwoStepReconstructionParameters
Base.@kwdef struct SinglePatchTwoStepReconstructionParameters{L_H, L_L<:AbstractSystemMatrixLoadingParameter, S<:AbstractLinearSolver, SP_H, SP_L<:AbstractSolverIterationParameters, R_H, R_L<:AbstractRegularization} <: AbstractSinglePatchReconstructionParameters
  # Threshhold
  Γ::Float64
  # File
  sf::MPIFile
  sfLoadHigh::L_H
  sfLoadLow::L_L
  # Solver
  solver::Type{S}
  solverParams_high::SP_H
  solverParams_low::SP_L
  reg_high::Vector{R_H} = AbstractRegularization[]
  reg_low::Vector{R_L} = AbstractRegularization[]
end
Base.@kwdef mutable struct SinglePatchTwoStepReconstructionAlgorithm{PR, R, PT} <: AbstractMPIReconstructionAlgorithm where {PR, R, PT<:AbstractMPIRecoParameters}
  params::SinglePatchParameters{PR, R, PT}
  # Could also do reconstruction progress meter here
  algoHigh::SinglePatchReconstructionAlgorithm
  algoLow::SinglePatchReconstructionAlgorithm
  output::Channel{Any}
end

# Bit hacky: Create transparent parameter to give to inner algorithm
Base.@kwdef mutable struct TwoStepSubstractionPreProcessingParameter{PR<:AbstractPreProcessingParameters} <: AbstractPreProcessingParameters
  pre::PR
  proj::Vector{ComplexF32} = zeros(ComplexF32, 1)
end
function Base.getproperty(param::TwoStepSubstractionPreProcessingParameter, field::Symbol)
  if field == :proj || field == :pre
    return getfield(param, field)
  else
    return getfield(getfield(param, :pre), field)
  end
end

function RecoUtils.process(t::Type{<:AbstractMPIReconstructionAlgorithm}, f::MPIFile, params::TwoStepSubstractionPreProcessingParameter)
  meas = process(t, f, params.pre)
  return meas - params.proj
end

function SinglePatchReconstruction(params::SinglePatchParameters{<:CommonPreProcessingParameters, <:SinglePatchTwoStepReconstructionParameters, PT}) where {PT <:AbstractPostProcessingParameters}
  return SinglePatchTwoStepReconstructionAlgorithm(params)
end
function SinglePatchTwoStepReconstructionAlgorithm(params::SinglePatchParameters{<:CommonPreProcessingParameters, <:SinglePatchTwoStepReconstructionParameters, PT}) where {PT <:AbstractPostProcessingParameters}
  recoHigh = SinglePatchReconstructionParameter(; sf = params.reco.sf, sfLoad = params.reco.sfLoadHigh, solver = params.reco.solver, solverParams = params.reco.solverParams_high, reg = params.reco.reg_high)
  recoLow = SinglePatchReconstructionParameter(; sf = params.reco.sf, sfLoad = params.reco.sfLoadLow, solver = params.reco.solver, solverParams = params.reco.solverParams_low, reg = params.reco.reg_low)
  algoHigh = SinglePatchReconstruction(SinglePatchParameters(params.pre, recoHigh, params.post))
  # First load proper S, grid and freqs
  algoLow = SinglePatchReconstruction(SinglePatchParameters(params.pre, recoLow, params.post))
  # Then construct "custom" SinglePatchAlgorithm
  paramsLow = SinglePatchParameters(TwoStepSubstractionPreProcessingParameter(;pre = algoLow.params.pre), algoLow.params.reco, algoLow.params.post)
  algoLow = SinglePatchReconstructionAlgorithm(paramsLow, params.reco.sf, algoLow.S, algoLow.grid, algoLow.freqs, algoLow.output)
  return SinglePatchTwoStepReconstructionAlgorithm(params, algoHigh, algoLow, Channel{Any}(Inf))
end

RecoUtils.take!(algo::SinglePatchTwoStepReconstructionAlgorithm) = Base.take!(algo.output)

function RecoUtils.put!(algo::SinglePatchTwoStepReconstructionAlgorithm, data::MPIFile)
  # First reco
  cPre = reconstruct(algo.algoHigh, data)

  # Thresholding
  thresh = maximum(abs.(cPre))*algo.params.reco.Γ
  cThresh = map(x-> abs(x) < thresh ? 0.0 : x, cPre.data)

  # Projection into raw data space
  uProj = map(ComplexF32,algo.algoLow.S*vec(cThresh))

  # Prepare subtraction
  algo.algoLow.params.pre.proj = uProj

  # Second reconstruction
  cPost = reconstruct(algo.algoLow, data)

  # Addition
  result = cPost + cThresh

  Base.put!(algo.output, result)
end

function RecoUtils.similar(algo::SinglePatchTwoStepReconstructionAlgorithm, data::MPIFile)
  algoHigh = RecoUtils.similar(algo.algoHigh, data)
  pre = fromKwargs(CommonPreProcessingParameters; toKwargs(algoHigh.params.pre)...)
  reco = RecoUtils.similar(algo, data, algo.params.reco)
  post = RecoUtils.similar(algo, data, algo.params.post)
  return SinglePatchReconstruction(SinglePatchParameters(pre, reco, post))
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
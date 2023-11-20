export SinglePatchTwoStepReconstructionParameters, SinglePatchTwoStepReconstructionAlgorithm
Base.@kwdef struct SinglePatchTwoStepReconstructionParameters{L_H<:AbstractSystemMatrixLoadingParameter, L_L<:AbstractSystemMatrixLoadingParameter, S<:AbstractLinearSolver, SP_H<:AbstractSolverParameters, SP_L<:AbstractSolverParameters, R_H<:AbstractRegularization, R_L<:AbstractRegularization} <: AbstractSinglePatchReconstructionParameters
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
Base.@kwdef mutable struct SinglePatchTwoStepReconstructionAlgorithm{P} <: AbstractSinglePatchReconstructionAlgorithm where {P<:AbstractSinglePatchAlgorithmParameters}
  params::P
  # Could also do reconstruction progress meter here
  algoHigh::SinglePatchReconstructionAlgorithm
  algoLow::SinglePatchReconstructionAlgorithm
  output::Channel{Any}
end

# Bit hacky: Create transparent parameter to give to inner algorithm
Base.@kwdef mutable struct TwoStepSubstractionPreProcessingParameter{B, PR<:AbstractMPIPreProcessingParameters{B}} <: AbstractMPIPreProcessingParameters{B}
  pre::PR
  proj::Array{ComplexF32} = zeros(ComplexF32, 0, 0)
end
function Base.getproperty(param::TwoStepSubstractionPreProcessingParameter, field::Symbol)
  if field == :proj || field == :pre
    return getfield(param, field)
  else
    return getproperty(getfield(param, :pre), field)
  end
end
function AbstractImageReconstruction.toKwargs(param::TwoStepSubstractionPreProcessingParameter; kwargs...)
  result = toKwargs(param.pre; kwargs...)
  result[:proj] = param.proj
  return result
end
function AbstractImageReconstruction.fromKwargs(type::Type{TwoStepSubstractionPreProcessingParameter{B, T}}; kwargs...) where {B, T}
  args = Dict{Symbol, Any}()
  dict = values(kwargs)
  for field in fieldnames(T)
    if haskey(dict, field)
      args[field] = getproperty(dict, field)
    end
  end
  proj = getproperty(dict, :proj)
  inner = T(;args...)
  return type(;pre = inner, proj = proj)
end

function process(t::Type{<:AbstractMPIRecoAlgorithm}, params::TwoStepSubstractionPreProcessingParameter, args...)
  meas = process(t, params.pre, args...)
  return meas .- params.proj
end

# Specify multi threaded variant, s.t. substraction happens on collected results
function process(algo::T, params::TwoStepSubstractionPreProcessingParameter, threadedInput::MultiThreadedInput, frequencies::Vector{CartesianIndex{2}}) where T <: SinglePatchReconstructionAlgorithm
  meas = process(algo, params.pre, threadedInput, frequencies)
  return meas .- params.proj
end

function SinglePatchReconstruction(params::SinglePatchParameters{<:CommonPreProcessingParameters, <:SinglePatchTwoStepReconstructionParameters, PT}) where {PT <:AbstractMPIPostProcessingParameters}
  return SinglePatchTwoStepReconstructionAlgorithm(params)
end
function SinglePatchTwoStepReconstructionAlgorithm(params::SinglePatchParameters{<:CommonPreProcessingParameters, <:SinglePatchTwoStepReconstructionParameters, PT}) where {PT <:AbstractMPIPostProcessingParameters}
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
AbstractImageReconstruction.parameter(algo::SinglePatchTwoStepReconstructionAlgorithm) = algo.params
AbstractImageReconstruction.take!(algo::SinglePatchTwoStepReconstructionAlgorithm) = Base.take!(algo.output)

function AbstractImageReconstruction.put!(algo::SinglePatchTwoStepReconstructionAlgorithm, data)
  # First reco
  cPre = reconstruct(algo.algoHigh, data)

  # Thresholding
  thresh = maximum(abs.(cPre))*algo.params.reco.Γ
  cThresh = map(x-> abs(x) < thresh ? 0.0 : x, cPre.data)

  # Projection into raw data space
  uProj = zeros(ComplexF32, size(algo.algoLow.S, 1), 1, size(cThresh)[end])
  cMapped = map(ComplexF32, cThresh)
  for frame = 1:size(uProj, 3)
    uProj[:, 1, frame] = algo.algoLow.S*vec(cMapped[1, :, :, :, frame])
  end

  # Prepare subtraction
  algo.algoLow.params.pre.proj = uProj

  # Second reconstruction
  cPost = reconstruct(algo.algoLow, data)

  # Addition
  result = AxisArray(cPost.data + cThresh, cPost.data.axes)
  result = ImageMeta(result, properties(cPost))

  Base.put!(algo.output, result)
end
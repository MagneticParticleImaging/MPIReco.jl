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


function SinglePatchReconstruction(params::SinglePatchParameters{<:CommonPreProcessingParameters, <:SinglePatchTwoStepReconstructionParameters, PT}) where {PT <:AbstractPostProcessingParameters}
  recoHigh = SinglePatchReconstructionParameter(; sf = params.reco.sf, sfLoad = params.reco.sfLoadHigh, solver = params.reco.solver, solverParams = params.reco.solverParams_high, reg = params.reco.reg_high)
  recoLow = SinglePatchReconstructionParameter(; sf = params.reco.sf, sfLoad = params.reco.sfLoadLow, solver = params.reco.solver, solverParams = params.reco.solverParams_low, reg = params.reco.reg_low)
  algoHigh = SinglePatchReconstruction(SinglePatchParameters(params.pre, recoHigh, params.post))
  algoLow = SinglePatchReconstruction(SinglePatchParameters(params.pre, recoLow, params.post))
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

  # Subtraction
  uMeas_low = process(algo.algoLow, data, algo.algoLow.params.pre)
  uCorr = uMeas_low - uProj

  # Second reconstruction
  temp = process(algo.algoLow, uCorr, algo.algoLow.params.reco)
  temp = process(algo.algoLow, temp, algo.algoLow.params.post)
  pixspacing = (spacing(algo.algoLow.grid) ./ acqGradient(data)[1] .* acqGradient(algo.algoLow.params.reco.sf)[1])*1000u"mm"
  offset = (ffPos(data) .- 0.5 .* calibFov(algo.algoLow.params.reco.sf))*1000u"mm" .+ 0.5 .* pixspacing
  dt = acqNumAverages(data)*dfCycle(data)*algo.algoLow.params.pre.numAverages*1u"s"
  axis = ImageAxisParameter(pixspacing=pixspacing, offset=offset, dt = dt)
  imMeta = ImageMetadataSystemMatrixParameter(data, algo.algoLow.params.reco.sf, algo.algoLow.grid, axis)
  cPost = process(AbstractMPIReconstructionAlgorithm, temp, imMeta)

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
Base.@kwdef struct SinglePatchReconstructionParameter{L<:AbstractSystemMatrixLoadingParameter, S<:AbstractLinearSolver,
   SP<:AbstractSolverParameters, R<:AbstractRegularization, W<:AbstractWeightingParameters} <: AbstractSinglePatchReconstructionParameters
  # File
  sf::MPIFile
  sfLoad::L
  # Solver
  solver::Type{S}
  solverParams::SP
  reg::Vector{R} = AbstractRegularization[]
  weightingParams::W = NoWeightingParamters()
end

Base.@kwdef mutable struct SinglePatchReconstructionAlgorithm{P} <: AbstractSinglePatchReconstructionAlgorithm where {P<:AbstractSinglePatchAlgorithmParameters}
  params::P
  # Could also do reconstruction progress meter here
  sf::Union{MPIFile, Vector{MPIFile}}
  S::AbstractArray
  grid::RegularGridPositions
  freqs::Vector{CartesianIndex{2}}
  output::Channel{Any}
end

function SinglePatchReconstruction(params::SinglePatchParameters{<:AbstractMPIPreProcessingParameters, R, PT}) where {R<:AbstractSinglePatchReconstructionParameters, PT <:AbstractMPIPostProcessingParameters}
  return SinglePatchReconstructionAlgorithm(params)
end
function SinglePatchReconstructionAlgorithm(params::SinglePatchParameters{<:AbstractMPIPreProcessingParameters, R, PT}) where {R<:AbstractSinglePatchReconstructionParameters, PT <:AbstractMPIPostProcessingParameters}
  freqs, S, grid = prepareSystemMatrix(params.reco)
  return SinglePatchReconstructionAlgorithm(params, params.reco.sf, S, grid, freqs, Channel{Any}(Inf))
end
recoAlgorithmTypes(::Type{SinglePatchReconstruction}) = SystemMatrixBasedAlgorithm()
AbstractImageReconstruction.parameter(algo::SinglePatchReconstructionAlgorithm) = algo.origParam

function prepareSystemMatrix(reco::SinglePatchReconstructionParameter{L,S}) where {L<:AbstractSystemMatrixLoadingParameter, S<:AbstractLinearSolver}
  freqs, sf, grid = process(AbstractMPIRecoAlgorithm, reco.sfLoad, reco.sf)
  sf, grid = prepareSF(S, sf, grid) 
  return freqs, sf, grid
end


AbstractImageReconstruction.take!(algo::SinglePatchReconstructionAlgorithm) = Base.take!(algo.output)

function AbstractImageReconstruction.put!(algo::SinglePatchReconstructionAlgorithm, data::MPIFile)
  consistenceCheck(algo.sf, data)
  
  result = process(algo, algo.params, data, algo.freqs)
  
  # Create Image (maybe image parameter as post params?)
  # TODO make more generic to apply to other pre/reco params as well (pre.numAverage main issue atm)
  pixspacing = (spacing(algo.grid) ./ acqGradient(data)[1] .* acqGradient(algo.sf)[1])*1000u"mm"
  offset = (ffPos(data) .- 0.5 .* calibFov(algo.sf))*1000u"mm" .+ 0.5 .* pixspacing
  dt = acqNumAverages(data)*dfCycle(data)*algo.params.pre.numAverages*1u"s"
  im = makeAxisArray(result, pixspacing, offset, dt)
  result = ImageMeta(im, generateHeaderDict(algo.sf, data))

  Base.put!(algo.output, result)
end

function process(algo::SinglePatchReconstructionAlgorithm, params::AbstractMPIPreProcessingParameters, f::MPIFile, args...)
  result = process(typeof(algo), params, f, args...)
  if eltype(algo.S) != eltype(result)
    @warn "System matrix and measurement have different element data type. Mapping measurment data to system matrix element type."
    result = map(eltype(algo.S),result)
  end
  return result
end


function process(algo::SinglePatchReconstructionAlgorithm, params::SinglePatchReconstructionParameter, u::Array)
  weights = process(algo, params.weightingParams, u)

  B = getLinearOperator(algo, params)

  solver = LeastSquaresParameters(solver = params.solver, op = B, S = algo.S, reg = params.reg, solverParams = params.solverParams, weights = weights)

  result = process(algo, solver, u)

  return gridresult(result, algo.grid, algo.sf)
end

process(algo::SinglePatchReconstructionAlgorithm, params::ChannelWeightingParameters, u::Array) = map(real(eltype(algo.S)), process(typeof(algo), params, algo.freqs))

function getLinearOperator(algo::SinglePatchReconstructionAlgorithm, params::SinglePatchReconstructionParameter{<:DenseSystemMatixLoadingParameter, S}) where {S}
  return nothing
end

function getLinearOperator(algo::SinglePatchReconstructionAlgorithm, params::SinglePatchReconstructionParameter{<:SparseSystemMatrixLoadingParameter, S}) where {S}
  return createLinearOperator(params.sfLoad.sparseTrafo, eltype(algo.S); shape=tuple(shape(algo.grid)...))
end

### Multi-Threading
function AbstractImageReconstruction.put!(algo::SinglePatchReconstructionAlgorithm, input::MultiThreadedInput)
  data = input.inputs[1]
  consistenceCheck(algo.sf, data)
  
  result = process(algo, algo.params, input.scheduler, data)
  
  # Create Image (maybe image parameter as post params?)
  # TODO make more generic to apply to other pre/reco params as well (pre.numAverage main issue atm)
  pixspacing = (spacing(algo.grid) ./ acqGradient(data)[1] .* acqGradient(algo.sf)[1])*1000u"mm"
  offset = (ffPos(data) .- 0.5 .* calibFov(algo.sf))*1000u"mm" .+ 0.5 .* pixspacing
  dt = acqNumAverages(data)*dfCycle(data)*algo.params.pre.numAverages*1u"s"
  im = makeAxisArray(result, pixspacing, offset, dt)
  result = ImageMeta(im, generateHeaderDict(algo.sf, data))

  Base.put!(algo.output, result)
end

function process(algo::SinglePatchReconstructionAlgorithm, params::SinglePatchParameters, scheduler::AbstractMultiThreadedProcessing, data::MPIFile)
  frames = collect(params.pre.pre.frames) 
  # TODO consider averaging and periodgrouping
  for frame in frames
    pre = fromKwargs(CommonPreProcessingParameters; toKwargs(params.pre)..., frames = [frame])
    p = SinglePatchParameters(pre, reco = params.reco, post = params.post)
    put!(scheduler, algo, p, data)
  end
  result = nothing
  for frame in frames
    if isnothing(result)
      result = take!(scheduler)
    else
      result = cat(result, take!(scheduler); dims = 5)
    end
  end
  return result
end
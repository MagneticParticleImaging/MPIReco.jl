Base.@kwdef struct SinglePatchReconstructionParameter{L<:AbstractSystemMatrixLoadingParameter, S<:AbstractLinearSolver,
   SP<:AbstractSolverParameters, R<:AbstractRegularization} <: AbstractSinglePatchReconstructionParameters
  # File
  sf::MPIFile
  sfLoad::L
  # Solver
  solver::Type{S}
  solverParams::SP
  reg::Vector{R} = AbstractRegularization[]
  # weightingType::WeightingType = WeightingType.None
end

Base.@kwdef mutable struct SinglePatchReconstructionAlgorithm{P} <: AbstractSinglePatchReconstructionAlgorithm where {P<:AbstractSinglePatchAlgorithmParameters}
  params::P
  # Could also do reconstruction progress meter here
  origParam::Union{AbstractSinglePatchAlgorithmParameters, Nothing} = nothing
  sf::Union{MPIFile, Vector{MPIFile}}
  S::AbstractArray
  grid::RegularGridPositions
  freqs::Vector{Int64}
  output::Channel{Any}
end

function SinglePatchReconstruction(params::SinglePatchParameters{<:AbstractMPIPreProcessingParameters, R, PT}) where {R<:AbstractSinglePatchReconstructionParameters, PT <:AbstractMPIPostProcessingParameters}
  return SinglePatchReconstructionAlgorithm(params)
end
function SinglePatchReconstructionAlgorithm(params::SinglePatchParameters{<:AbstractMPIPreProcessingParameters, R, PT}) where {R<:AbstractSinglePatchReconstructionParameters, PT <:AbstractMPIPostProcessingParameters}
  freqs, S, grid = prepareSystemMatrix(params.reco)
  filter = FrequencyFilteredPreProcessingParameters(freqs, params.pre)
  filteredParams = SinglePatchParameters(filter, params.reco, params.post)
  return SinglePatchReconstructionAlgorithm(filteredParams, params, params.reco.sf, S, grid, freqs, Channel{Any}(Inf))
end
recoAlgorithmTypes(::Type{SinglePatchReconstruction}) = SystemMatrixBasedAlgorithm()
AbstractImageReconstruction.parameter(algo::SinglePatchReconstructionAlgorithm) = algo.origParam

function prepareSystemMatrix(reco::SinglePatchReconstructionParameter{L,S}) where {L<:AbstractSystemMatrixLoadingParameter, S<:AbstractLinearSolver}
  freqs, sf, grid = process(AbstractMPIRecoAlgorithm, reco.sf, reco.sfLoad)
  sf, grid = prepareSF(S, sf, grid) 
  return freqs, sf, grid
end


AbstractImageReconstruction.take!(algo::SinglePatchReconstructionAlgorithm) = Base.take!(algo.output)

function AbstractImageReconstruction.put!(algo::SinglePatchReconstructionAlgorithm, data::MPIFile)
  consistenceCheck(algo.sf, data)
  
  result = process(algo, data, algo.params)
  
  # Create Image (maybe image parameter as post params?)
  # TODO make more generic to apply to other pre/reco params as well (pre.numAverage main issue atm)
  pixspacing = (spacing(algo.grid) ./ acqGradient(data)[1] .* acqGradient(algo.sf)[1])*1000u"mm"
  offset = (ffPos(data) .- 0.5 .* calibFov(algo.sf))*1000u"mm" .+ 0.5 .* pixspacing
  dt = acqNumAverages(data)*dfCycle(data)*algo.params.pre.numAverages*1u"s"
  im = makeAxisArray(result, pixspacing, offset, dt)
  result = ImageMeta(im, generateHeaderDict(algo.sf, data))

  Base.put!(algo.output, result)
end

function process(algo::SinglePatchReconstructionAlgorithm, f::MPIFile, params::AbstractMPIPreProcessingParameters)
  result = process(typeof(algo), f, params)
  if eltype(algo.S) != eltype(result)
    @warn "System matrix and measurement have different element data type. Mapping measurment data to system matrix element type."
    result = map(eltype(algo.S),result)
  end
  return result
end


function process(algo::SinglePatchReconstructionAlgorithm, u::Array, params::SinglePatchReconstructionParameter)
  weights = nothing # getWeights(...)

  B = getLinearOperator(algo, params)

  solver = LeastSquaresParameters(params.solver, B, algo.S, params.reg, params.solverParams)

  result = process(algo, u, solver)

  return gridresult(result, algo.grid, algo.sf)
end

function getLinearOperator(algo::SinglePatchReconstructionAlgorithm, params::SinglePatchReconstructionParameter{<:DenseSystemMatixLoadingParameter, S}) where {S}
  return nothing
end

function getLinearOperator(algo::SinglePatchReconstructionAlgorithm, params::SinglePatchReconstructionParameter{<:SparseSystemMatrixLoadingParameter, S}) where {S}
  return createLinearOperator(params.sfLoad.sparseTrafo, eltype(algo.S); shape=tuple(shape(algo.grid)...))
end
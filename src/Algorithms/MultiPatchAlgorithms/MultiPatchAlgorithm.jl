Base.@kwdef struct MultiPatchReconstructionParameter{L<:AbstractSystemMatrixLoadingParameter, S<:AbstractLinearSolver,
    SP<:AbstractSolverParameters, R<:AbstractRegularization} <: AbstractMultiPatchReconstructionParameters
   # File
   sf::MultiMPIFile
   sfLoad::L
   solverParams::SP
   reg::Vector{R} = AbstractRegularization[]
   # weightingType::WeightingType = WeightingType.None
 end
 
 Base.@kwdef mutable struct MultiPatchReconstructionAlgorithm{P} <: AbstractMultiPatchReconstructionAlgorithm where {P<:AbstractMultiPatchAlgorithmParameters}
   params::P
   # Could also do reconstruction progress meter here
   origParam::Union{AbstractMultiPatchAlgorithmParameters, Nothing} = nothing
   sf::MultiMPIFile
   S::AbstractArray
   grid::RegularGridPositions
   freqs::Vector{Int64}
   output::Channel{Any}
 end
 
 function MultiPatchReconstruction(params::MultiPatchParameters{<:AbstractPreProcessingParameters, R, PT}) where {R<:AbstractMultiPatchReconstructionParameters, PT <:AbstractPostProcessingParameters}
   return MultiPatchReconstructionAlgorithm(params)
 end
 function MultiPatchReconstructionAlgorithm(params::MultiPatchParameters{<:AbstractPreProcessingParameters, R, PT}) where {R<:AbstractMultiPatchReconstructionParameters, PT <:AbstractPostProcessingParameters}
   freqs, S, grid = prepareSystemMatrix(params.reco)
   filter = fromKwargs(FrequencyFilteredPreProcessingParameters; frequencies = freqs, toKwargs(params.pre; flatten = DataType[])...)
   filteredParams = MultiPatchParameters(filter, params.reco, params.post)
   return MultiPatchReconstructionAlgorithm(filteredParams, params, params.reco.sf, S, grid, freqs, Channel{Any}(Inf))
 end
 recoAlgorithmTypes(::Type{MultiPatchReconstruction}) = SystemMatrixBasedAlgorithm()
 RecoUtils.parameter(algo::MultiPatchReconstructionAlgorithm) = algo.origParam
 
 function prepareSystemMatrix(reco::MultiPatchReconstructionParameter{L,S}) where {L<:AbstractSystemMatrixLoadingParameter, S<:AbstractLinearSolver}
   freqs, sf, grid = process(AbstractMPIReconstructionAlgorithm, reco.sf, reco.sfLoad)
   sf, grid = prepareSF(S, sf, grid) 
   return freqs, sf, grid
 end
  
 RecoUtils.take!(algo::MultiPatchReconstructionAlgorithm) = Base.take!(algo.output)
 
 function RecoUtils.put!(algo::MultiPatchReconstructionAlgorithm, data::MPIFile)
   consistenceCheck(algo.sf, data)

   algo.FFOp = process(algo, data, algo.opParams)
   
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
 
 function RecoUtils.process(algo::MultiPatchReconstructionAlgorithm, f::MPIFile, params::AbstractPreProcessingParameters)
   result = process(typeof(algo), f, params)
   if eltype(algo.S) != eltype(result)
     @warn "System matrix and measurement have different element data type. Mapping measurment data to system matrix element type."
     result = map(eltype(algo.S),result)
   end
   return result
 end
 
 
 function RecoUtils.process(algo::MultiPatchReconstructionAlgorithm, u::Array, params::MultiPatchReconstructionParameter)
   solver = LeastSquaresParameters(Kaczmarz, nothing, algo.FFOp, params.reg, params.solverParams)
 
   result = process(algo, u, solver)
 
   return gridresult(result, algo.grid, algo.sf)
 end
export SinglePatchTemporalRegularizationAlgorithm, SinglePatchTemporalRegularizationReconstructionParameter
Base.@kwdef struct SinglePatchTemporalRegularizationReconstructionParameter{L<:DenseSystemMatixLoadingParameter,
  SP<:AbstractSolverParameters} <: AbstractSinglePatchReconstructionParameters
  # File
  sf::MPIFile
  sfLoad::L
  # Solver
  solverParams::SP
  λ::Float32
  β::Float32
  # weightingType::WeightingType = WeightingType.None
  idxFG::Union{Nothing, UnitRange{Int64}, Vector{Int64}} = nothing
  idxFG::Union{Nothing, UnitRange{Int64}, Vector{Int64}} = nothing
  bgDict::BGDictParameter
end

Base.@kwdef mutable struct SinglePatchTemporalRegularizationAlgorithm{P} <: AbstractSinglePatchReconstructionAlgorithm where {P<:AbstractSinglePatchAlgorithmParameters}
  params::P
  origParam::Union{AbstractSinglePatchAlgorithmParameters,Nothing} = nothing
  sf::Union{MPIFile,Vector{MPIFile}}
  S::AbstractArray
  bgDict::AbstractArray
  idxFG::Union{Nothing, UnitRange{Int64}, Vector{Int64}} = nothing
  idxFG::Union{Nothing, UnitRange{Int64}, Vector{Int64}} = nothing
  grid::RegularGridPositions
  freqs::Vector{Int64}
  output::Channel{Any}
end

function SinglePatchReconstruction(params::SinglePatchParameters{<:AbstractMPIPreProcessingParameters,<:SinglePatchTemporalRegularizationReconstructionParameter,PT}) where {PT<:AbstractMPIPostProcessingParameters}
  return SinglePatchTemporalRegularizationAlgorithm(params)
end
function SinglePatchTemporalRegularizationAlgorithm(params::SinglePatchParameters{<:AbstractMPIPreProcessingParameters,R,PT}) where {R<:SinglePatchTemporalRegularizationReconstructionParameter,PT<:AbstractMPIPostProcessingParameters}
  freqs, S, grid = prepareSystemMatrix(params.reco)
  filter = FrequencyFilteredPreProcessingParameters(freqs, params.pre)
  filteredParams = SinglePatchParameters(filter, params.reco, params.post)
  return SinglePatchTemporalRegularizationAlgorithm(filteredParams, params, params.reco.sf, S, process(SinglePatchTemporalRegularizationAlgorithm, freqs, params.reco.bgDict)
    ,params.reco.idxFG, params.reco.idxBG, grid, freqs, Channel{Any}(Inf))
end
recoAlgorithmTypes(::Type{SinglePatchTemporalRegularizationAlgorithm}) = SystemMatrixBasedAlgorithm()
AbstractImageReconstruction.parameter(algo::SinglePatchTemporalRegularizationAlgorithm) = algo.origParam

function prepareSystemMatrix(reco::SinglePatchTemporalRegularizationReconstructionParameter{L}) where {L<:AbstractSystemMatrixLoadingParameter}
  freqs, sf, grid = process(AbstractMPIRecoAlgorithm, reco.sf, reco.sfLoad)
  sf, grid = prepareSF(Kaczmarz, sf, grid)
  return freqs, sf, grid
end

AbstractImageReconstruction.take!(algo::SinglePatchTemporalRegularizationAlgorithm) = Base.take!(algo.output)

function AbstractImageReconstruction.put!(algo::SinglePatchTemporalRegularizationAlgorithm, data::MPIFile)
  consistenceCheck(algo.sf, data)

  result = process(algo, data, algo.params)

  # Create Image (maybe image parameter as post params?)
  # TODO make more generic to apply to other pre/reco params as well (pre.numAverage main issue atm)
  pixspacing = (spacing(algo.grid) ./ acqGradient(data)[1] .* acqGradient(algo.sf)[1]) * 1000u"mm"
  offset = (ffPos(data) .- 0.5 .* calibFov(algo.sf)) * 1000u"mm" .+ 0.5 .* pixspacing
  dt = acqNumAverages(data) * dfCycle(data) * algo.params.pre.numAverages * 1u"s"
  im = makeAxisArray(result, pixspacing, offset, dt)
  result = ImageMeta(im, generateHeaderDict(algo.sf, data))

  Base.put!(algo.output, result)
end


function process(algo::SinglePatchTemporalRegularizationAlgorithm, f::MPIFile, params::AbstractMPIPreProcessingParameters)
  result = process(typeof(algo), f, params)
  if eltype(algo.S) != eltype(result)
    @warn "System matrix and measurement have different element data type. Mapping measurment data to system matrix element type."
    result = map(eltype(algo.S), result)
  end
  return result
end


function process(algo::SinglePatchTemporalRegularizationAlgorithm, u::Array, params::SinglePatchTemporalRegularizationReconstructionParameter)
  weights = nothing # getWeights(...)

  L = size(u)[end]
  idxFG = isnothing(algo.idxFG) ? (collect(1:L)) : algo.idxFG
  idxFG = isnothing(algo.idxBG) ? (collect(1:L)) : algo.idxBG

  # Prepare Regularization
  reg = L2Regularization(Float32(params.λ))
  λ = RegularizedLeastSquares.normalize(params.solverParams.normalizeReg, reg, algo.S, nothing) # Scaling happens outside
  reg = L2Regularization(Float32(1))

  # TODO interpMeth as parameter
  op = TemporalRegularizationOperator(Float32(1/(sqrt(λ)))*S, Float32(1/(sqrt(β)))*bgDict, L; :linear, idxFG, idxBG )

  MSub, NSub, Q, J, Γ = getSizes(op)

  constraintMask = zeros(Bool, op.N)
  constraintMask[1:(NSub*J)] .= 1
  
  # Enforce no normalization
  solverParams = fromKwargs(typeof(params.solverParams); toKwargs(params.solverParams, overwrite = Dict{Symbol, Any}(:normalizeReg => NoNormalization()))...)
  solverParams = ConstraintMaskedSolverParameters(;constraintMask = constraintMask, params = params.solverParams)

  solver = LeastSquaresParameters(Kaczmarz, nothing, op, [reg], solverParams)

  temp = process(algo, u, solver)
  temp = real.( reshape(temp[1:(NSub*J),:],NSub,J) ./ sqrt(λ) ) 

  cInterp = similar(temp, size(c,1), op.L)
  for l=1:op.L
    for κ=1:length(op.idxCoeffsFG[l])
      cInterp[:,l] .+= op.coeffsFG[l][κ] * temp[:, op.idxCoeffsFG[l][κ]]
    end
  end

  return gridresult(result, algo.grid, algo.sf)
end
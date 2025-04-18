export SinglePatchTemporalRegularizationAlgorithm, SinglePatchTemporalRegularizationReconstructionParameter
Base.@kwdef struct SinglePatchTemporalRegularizationReconstructionParameter{L<:DenseSystemMatixLoadingParameter,
  SP<:AbstractSolverParameters, arrT <: AbstractArray} <: AbstractSinglePatchReconstructionParameters
  # File
  sf::MPIFile
  sfLoad::Union{L, AbstractUtilityReconstructionParameters{L}}
  arrayType::Type{arrT} = Array
  # Solver
  solverParams::SP
  λ::Float32
  β::Float32
  # weightingType::WeightingType = WeightingType.None
  idxFG::Union{Nothing, UnitRange{Int64}, Vector{Int64}} = nothing
  idxFG::Union{Nothing, UnitRange{Int64}, Vector{Int64}} = nothing
  bgDict::BGDictParameter
end

Base.@kwdef mutable struct SinglePatchTemporalRegularizationAlgorithm{P, arrT <: AbstractArray} <: AbstractSinglePatchReconstructionAlgorithm where {P<:AbstractSinglePatchAlgorithmParameters}
  params::P
  sf::Union{MPIFile,Vector{MPIFile}}
  S::AbstractArray
  arrayType::Type{arrT}
  bgDict::AbstractArray
  idxFG::Union{Nothing, UnitRange{Int64}, Vector{Int64}} = nothing
  idxFG::Union{Nothing, UnitRange{Int64}, Vector{Int64}} = nothing
  grid::RegularGridPositions
  freqs::Vector{CartesianIndex{2}}
  output::Channel{Any}
end

function SinglePatchReconstruction(params::SinglePatchParameters{<:AbstractMPIPreProcessingParameters,<:SinglePatchTemporalRegularizationReconstructionParameter,PT}) where {PT<:AbstractMPIPostProcessingParameters}
  return SinglePatchTemporalRegularizationAlgorithm(params)
end
function SinglePatchTemporalRegularizationAlgorithm(params::SinglePatchParameters{<:AbstractMPIPreProcessingParameters,R,PT}) where {R<:SinglePatchTemporalRegularizationReconstructionParameter,PT<:AbstractMPIPostProcessingParameters}
  freqs, S, grid, arrayType = prepareSystemMatrix(params.reco)
  return SinglePatchTemporalRegularizationAlgorithm(params, params.reco.sf, S, arrayType, process(SinglePatchTemporalRegularizationAlgorithm, params.reco.bgDict, freqs)
    ,params.reco.idxFG, params.reco.idxBG, grid, freqs, Channel{Any}(Inf))
end
recoAlgorithmTypes(::Type{SinglePatchTemporalRegularizationAlgorithm}) = SystemMatrixBasedAlgorithm()
AbstractImageReconstruction.parameter(algo::SinglePatchTemporalRegularizationAlgorithm) = algo.origParam

function prepareSystemMatrix(reco::SinglePatchTemporalRegularizationReconstructionParameter{L}) where {L<:AbstractSystemMatrixLoadingParameter}
  freqs, sf, grid = process(AbstractMPIRecoAlgorithm, reco.sfLoad, reco.sf, Kaczmarz, reco.arrayType)
  return freqs, sf, grid, reco.arrayType
end

Base.lock(algo::SinglePatchTemporalRegularizationAlgorithm) = lock(algo.output)
Base.unlock(algo::SinglePatchTemporalRegularizationAlgorithm) = unlock(algo.output)
Base.isready(algo::SinglePatchTemporalRegularizationAlgorithm) = isready(algo.output)
Base.wait(algo::SinglePatchTemporalRegularizationAlgorithm) = wait(algo.output)
AbstractImageReconstruction.take!(algo::SinglePatchTemporalRegularizationAlgorithm) = Base.take!(algo.output)

function process(algo::SinglePatchTemporalRegularizationAlgorithm, params::AbstractMPIPreProcessingParameters, f::MPIFile)
  result = process(typeof(algo), f, params)
  if eltype(algo.S) != eltype(result)
    @warn "System matrix and measurement have different element data type. Mapping measurment data to system matrix element type."
    result = map(eltype(algo.S), result)
  end
  result = adapt(algo.arrayType, result)
  return result
end


function process(algo::SinglePatchTemporalRegularizationAlgorithm, params::SinglePatchTemporalRegularizationReconstructionParameter, u::Array)
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

  solver = LeastSquaresParameters(solver = Kaczmarz, S = op, reg = [reg], solverParams = solverParams)

  temp = process(algo, solver, u)
  temp = real.( reshape(temp[1:(NSub*J),:],NSub,J) ./ sqrt(λ) ) 

  cInterp = similar(temp, size(c,1), op.L)
  for l=1:op.L
    for κ=1:length(op.idxCoeffsFG[l])
      cInterp[:,l] .+= op.coeffsFG[l][κ] * temp[:, op.idxCoeffsFG[l][κ]]
    end
  end

  return gridresult(result, algo.grid, algo.sf)
end
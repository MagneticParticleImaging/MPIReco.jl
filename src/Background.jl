export getBackgroundDictionaryComplete

export NoBackgroundCorrectionParameters
struct NoBackgroundCorrectionParameters <: AbstractMPIBackgroundCorrectionParameters end

export InternalBackgroundCorrectionParameters
Base.@kwdef struct InternalBackgroundCorrectionParameters <: AbstractMPIBackgroundCorrectionParameters
  interpolateBG::Bool = false
end

export ExternalBackgroundCorrection
abstract type ExternalBackgroundCorrection <: AbstractMPIBackgroundCorrectionParameters end

Base.@kwdef struct ExternalPreProcessedBackgroundCorrectionParameters{T} <: AbstractMPIBackgroundCorrectionParameters where {T<:ExternalBackgroundCorrection}
  numPeriodAverages::Int64 = 1
  numPeriodGrouping::Int64 = 1
  spectralLeakageCorrection::Bool = false
  loadasreal::Bool = false
  bgParams::T
end
# TODO Shorter struct names?

export SimpleExternalBackgroundCorrectionParameters
Base.@kwdef struct SimpleExternalBackgroundCorrectionParameters <: ExternalBackgroundCorrection
  emptyMeas::MPIFile
  bgFrames::Union{Vector{Int64}, UnitRange{Int64}} = [1]
end
function process(::Type{<:AbstractMPIRecoAlgorithm}, params::SimpleExternalBackgroundCorrectionParameters, data::Array)
  kwargs = toKwargs(params, overwrite = Dict{Symbol, Any}(:frames => params.bgFrames))
  empty = getMeasurementsFD(params.emptyMeas, false; bgCorrection = false, numAverages=length(bgFrames), kwargs...)
  return data .-empty
end
function process(::Type{<:AbstractMPIRecoAlgorithm}, params::ExternalPreProcessedBackgroundCorrectionParameters{SimpleExternalBackgroundCorrectionParameters}, data::Array, frequencies::Union{Vector{CartesianIndex{2}}, Nothing} = nothing)
  kwargs = toKwargs(params, overwrite = Dict{Symbol, Any}(:frames => params.bgParams.bgFrames), ignore = [:bgParams])
  empty = getMeasurementsFD(params.bgParams.emptyMeas, false; bgCorrection = false, numAverages=length(params.bgParams.bgFrames), kwargs..., frequencies = frequencies)
  return data .-empty
end

export LinearInterpolatedExternalBackgroundCorrectionParameters
Base.@kwdef struct LinearInterpolatedExternalBackgroundCorrectionParameters <: AbstractMPIBackgroundCorrectionParameters
  emptyMeas::MPIFile
  bgFrames::UnitRange{Int64} = 1:1
  bgFramesPost::UnitRange{Int64} = 1:1
end
function process(::Type{<:AbstractMPIRecoAlgorithm}, params::LinearInterpolatedExternalBackgroundCorrectionParameters, data::Array)
  kwargs = toKwargs(params)
  kwargs[:frames] = params.bgFrames
  empty = getMeasurementsFD(bgParams.emptyMeas, false; bgCorrection = false, numAverages=length(params.bgFrames), kwargs...)
  kwargs[:frames] = params.bgFramesPost
  emptyPost = getMeasurementsFD(params.emptyMeas, false; bgCorrection = false, numAverages=length(params.bgFramesPost), kwargs...)
  for l=1:size(result, 4)
    alpha = (l - mean(params.bgFrames)) / (mean(params.bgFramesPost) - mean(params.bgFrames))
    result[:,:,l] .-=  (1-alpha).*empty[:,:,1] .+ alpha.*emptyPost[:,:,1]
  end
  return result
end
function process(::Type{<:AbstractMPIRecoAlgorithm}, params::ExternalPreProcessedBackgroundCorrectionParameters{LinearInterpolatedExternalBackgroundCorrectionParameters}, data::Array, frequencies::Union{Vector{CartesianIndex{2}}, Nothing} = nothing)
  kwargs = toKwargs(params, ignore = [:bgParams])
  bgParams = params.bgParams
  kwargs[:frames] = bgParams.bgFrames
  empty = getMeasurementsFD(bgParams.emptyMeas, false; bgCorrection = false, numAverages=length(bgParams.bgFrames), kwargs..., frequencies = frequencies)
  kwargs[:frames] = bgParams.bgFramesPost
  emptyPost = getMeasurementsFD(bgParams.emptyMeas, false; bgCorrection = false, numAverages=length(bgParams.bgFramesPost), kwargs..., frequencies = frequencies)
  for l=1:size(result, 4)
    alpha = (l - mean(bgParams.bgFrames)) / (mean(bgParams.bgFramesPost) - mean(bgParams.bgFrames))
    result[:,:,l] .-=  (1-alpha).*empty[:,:,1] .+ alpha.*emptyPost[:,:,1]
  end
  return result
end

export AbstractBGDictLoader
abstract type AbstractBGDictLoader <: AbstractMPIRecoParameters end
export MeasurementBGDictLoader
Base.@kwdef struct MeasurementBGDictLoader{T} <: AbstractBGDictLoader where {T<:MPIFile}
  file::T
  bgFrames::Union{UnitRange{Int64}, Vector{Int64}} = measBGFrameIdx(file)
  numPeriodGrouping::Int64 = acqNumPeriodsPerFrame(file)
  numPeriodAverages::Int64 = 1
  bgAverages::Int64 = 1
end
function process(algoT::Type{<:AbstractMPIRecoAlgorithm}, params::MeasurementBGDictLoader, freqs::Vector{CartesianIndex{2}})
  uEmpty = getMeasurementsFD(params.file, false, frequencies=freqs, frames=params.bgFrames, numAverages=params.bgAverages, spectralLeakageCorrection=false, bgCorrection=false, numPeriodGrouping = params.numPeriodGrouping, numPeriodAverages = params.numPeriodAverages)
  return transpose(reshape(uEmpty, :, div(length(params.bgFrames),params.bgAverages)))
end
export SystemMatrixBGDictLoader
Base.@kwdef struct SystemMatrixBGDictLoader{T} <: AbstractBGDictLoader where {T<:MPIFile}
  file::T
end
function process(algoT::Type{<:AbstractMPIRecoAlgorithm}, params::SystemMatrixBGDictLoader, freqs::Vector{CartesianIndex{2}})
  idxBGFrames = measBGFrameIdx(params.file)
  D = measData(params.file, idxBGFrames)
  D_ = reshape(D, size(D,1), size(D,2), size(D,3), size(D,4))
  return reshape(D_[:,freqs,:],length(idxBGFrames),:)
end

function getBackgroundDictionaryComplete(fSF::MPIFile, f::MPIFile, frequencies,
                                       bgFrames=nothing, numBGAverages=1)
  idxBGFrames = measBGFrameIdx(fSF)

  D = measData(fSF, idxBGFrames)
  D_ = reshape(D, size(D,1), size(D,2)*size(D,3), size(D,4))
  bgdata = reshape(D_[:,frequencies,:],length(idxBGFrames),:)

  if bgFrames != nothing
    uEmpty = getMeasurementsFD(f, frequencies=frequencies, frames=bgFrames, numAverages=numBGAverages,
                               spectralLeakageCorrection=false, bgCorrection=false)
    bgdata2 = transpose(reshape(uEmpty, :, div(length(bgFrames),numBGAverages)))

    bgdata = cat(bgdata2,bgdata,dims=1)
  end

  return svd(transpose(bgdata))
end

export BGDictParameter
Base.@kwdef struct BGDictParameter <: AbstractMPIRecoParameters
  loader::Vector{AbstractBGDictLoader}
  dictSize::Int64 = 2
end

function process(algoT::Type{<:AbstractMPIRecoAlgorithm}, params::BGDictParameter, freqs)
  U,S,V = process(algoT, freqs, params.loader)
  for l=1:params.dictSize
    U[:,l] *= (S[l] / S[1])^(1/2)
  end
  return U[:,1:params.dictSize]
end

function process(algoT::Type{<:AbstractMPIRecoAlgorithm}, loader::Vector{T}, freqs) where {T<:AbstractBGDictLoader}
  tempBGs = []
  for load in loader
    push!(tempBGs, process(algoT, freqs, load))
  end
  return svd(transpose(cat(tempBGs..., dims=1)))
end

function getBackgroundDictionary(fSF::MPIFile, f::MPIFile, frequencies,
                                 bgDictSize::Int=2, bgFrames=nothing, numBGAverages=1)
  U,S,V = getBackgroundDictionaryComplete(fSF, f, frequencies, bgFrames, numBGAverages)

  for l=1:bgDictSize
    U[:,l] *= (S[l] / S[1])^(1/2)
  end
  #return transpose(bgdata)[:,1:bgDictSize]
  return U[:,1:bgDictSize]
end

getBackgroundDictionary(::Any, ::Any, ::Any, bgDictSize::Nothing, ::Any) = nothing



"""
Joint reconstruction of MPI signal and background. Implemented in a low-level
fashion
"""
function reconstruction(S, u::Array, bgDict::AbstractMatrix;
                        sparseTrafo = nothing, beta = 0.1, β=beta,
                        lambd=0.0, lambda=lambd, λ=lambda, progress=nothing,
                        solver = "Kaczmarz",
                        weights=nothing, enforceReal=false, enforcePositive=false,
                        relativeLambda=true, backgroundCoefficients = nothing, kargs...)

  N = size(S,2)
  M = div(length(S), N)

  L = size(u)[end]
  u = reshape(u, M, L)
  c = zeros(N,L)

  #β /= size(bgDict,2)#*1000

  if sum(abs.(λ)) > 0 && solver != "fusedlasso" && relativeLambda
    trace = calculateTraceOfNormalMatrix(S,weights)

    @debug "REL λ =  $(trace / N) "

    λ *= trace / N
    setlambda(S,λ)
  end

  G = transpose( cat(Float32(1/(sqrt(λ)))*transpose(S), Float32(1/(sqrt(β)))*transpose(bgDict), dims=1) )

  constraintMask = zeros(Bool, size(G,2))
  constraintMask[1:N] .= 1

  solv = createLinearSolver(solver, G; weights=weights, λ=1.0, constraintMask=constraintMask,
                            sparseTrafo=sparseTrafo, enforceReal=enforceReal,
			                enforcePositive=enforcePositive, kargs...)

  progress==nothing ? p = Progress(L, 1, "Reconstructing data...") : p = progress
  for l=1:L

    y = solve(solv, u[:,l])
    d = y[1:N,:] ./ sqrt(λ)

    if backgroundCoefficients != nothing
      append!(backgroundCoefficients, vec(y[(N+1):end,:] ./ sqrt(β)))
    end

    if sparseTrafo != nothing
      d[:] = sparseTrafo*d # backtrafo from dual space
    end

    c[:,l] = real( d )
    next!(p)
    sleep(0.001)
  end

  return c
end

export getBackgroundDictionaryComplete

export NoBackgroundCorrection
struct NoBackgroundCorrection <: AbstractBackgroundCorrectionParameters end

export InternalBackgroundCorrection
Base.@kwdef struct InternalBackgroundCorrection <: AbstractBackgroundCorrectionParameters
  interpolateBG::Bool = false
end

export ExternalBackgroundCorrection
abstract type ExternalBackgroundCorrection <: AbstractBackgroundCorrectionParameters end

Base.@kwdef struct ExternalPreProcessedBackgroundCorrectionParameters{T} <: AbstractBackgroundCorrectionParameters where {T<:ExternalBackgroundCorrection}
  numPeriodAverages::Int64 = 1
  numPeriodGrouping::Int64 = 1
  spectralLeakageCorrection::Bool = false
  loadasreal::Bool = false
  bgParams::T
end

Base.@kwdef struct FrequencyFilteredBackgroundCorrection{T} <: AbstractBackgroundCorrectionParameters where {T<:ExternalBackgroundCorrection}
  frequencies::Vector{Int64}
  numPeriodAverages::Int64 = 1
  numPeriodGrouping::Int64 = 1
  spectralLeakageCorrection::Bool = false
  loadasreal::Bool = false
  bgParams::T
end

export SimpleExternalBackgroundCorrection
Base.@kwdef struct SimpleExternalBackgroundCorrection <: AbstractBackgroundCorrectionParameters
  emptyMeas::MPIFile
  bgFrames::UnitRange{Int64} = 1:1
end
function RecoUtils.process(::Type{<:AbstractMPIReconstructionAlgorithm}, data::Array, params::SimpleExternalBackgroundCorrection)
  kwargs = toKwargs(params)
  kwargs[:frames] = params.bgFrames
  empty = getMeasurementsFD(bgParams.emptyMeas, false, bgCorrection = false, numAverages=length(bgFrames), kwargs...)
  return data .-empty
end
function RecoUtils.process(::Type{<:AbstractMPIReconstructionAlgorithm}, data::Array, fparams::FrequencyFilteredBackgroundCorrection{SimpleExternalBackgroundCorrection})
  params = fparams.bgParams
  kwargs = toKwargs(params)
  kwargs[:frames] = params.bgFrames
  empty = getMeasurementsFD(bgParams.emptyMeas, false, bgCorrection = false, numAverages=length(bgFrames), kwargs...)
  return data .-empty
end

export LinearInterpolatedExternalBackgroundCorrection
Base.@kwdef struct LinearInterpolatedExternalBackgroundCorrection <: AbstractBackgroundCorrectionParameters
  emptyMeas::MPIFile
  bgFrames::UnitRange{Int64} = 1:1
  bgFramesPost::UnitRange{Int64} = 1:1
end
function RecoUtils.process(::Type{<:AbstractMPIReconstructionAlgorithm}, data::Array, params::LinearInterpolatedExternalBackgroundCorrection)
  kwargs = toKwargs(params)
  kwargs[:frames] = params.bgFrames
  empty = getMeasurementsFD(bgParams.emptyMeas, false, bgCorrection = false, numAverages=length(params.bgFrames), kwargs...)
  kwargs[:frames] = params.bgFramesPost
  emptyPost = getMeasurementsFD(params.emptyMeas, false, bgCorrection = false, numAverages=length(params.bgFramesPost), kwargs...)
  for l=1:size(result, 4)
    alpha = (l - mean(params.bgFrames)) / (mean(params.bgFramesPost) - mean(params.bgFrames))
    result[:,:,l] .-=  (1-alpha).*empty[:,:,1] .+ alpha.*emptyPost[:,:,1]
  end
  return result
end
function RecoUtils.process(::Type{<:AbstractMPIReconstructionAlgorithm}, data::Array, fparams::FrequencyFilteredBackgroundCorrection{LinearInterpolatedExternalBackgroundCorrection})
  params = fparams.bgParams
  kwargs = toKwargs(params)
  kwargs[:frequencies] = fparams.freqs
  kwargs[:frames] = params.bgFrames
  empty = getMeasurementsFD(bgParams.emptyMeas, false, bgCorrection = false, numAverages=length(params.bgFrames), kwargs...)
  kwargs[:frames] = params.bgFramesPost
  emptyPost = getMeasurementsFD(params.emptyMeas, false, bgCorrection = false, numAverages=length(params.bgFramesPost), kwargs...)
  for l=1:size(result, 4)
    alpha = (l - mean(params.bgFrames)) / (mean(params.bgFramesPost) - mean(params.bgFrames))
    result[:,:,l] .-=  (1-alpha).*empty[:,:,1] .+ alpha.*emptyPost[:,:,1]
  end
  return result
end


export DictionaryBasedBackgroundCorrection
Base.@kwdef struct DictionaryBasedBackgroundCorrection <: AbstractBackgroundCorrectionParameters
  # TODO
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
                        solver = "kaczmarz",
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

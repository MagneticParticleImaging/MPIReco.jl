# Maybe two level struct: Params for introspection and then struct itself
mutable struct SimplePreProcessing <: AbstractMPIReconstructionAlgorithm
  # Params
  freq::Array
  frames
  numAverages
  loadasreal
  spectralLeakageCorrection
  bgCorrection
  numPeriodAverages
  numPeriodGrouping
  # Pipeline
  algo::AbstractMPIReconstructionAlgorithm
  outputChannel::Channel{Any}
end
function SimplePreProcessing(algo::AbstractMPIReconstructionAlgorithm; freq::Array, frames = nothing, 
  numAverages = 1,
  loadasreal = false,
  spectralLeakageCorrection = true,
  bgCorrection=false, numPeriodAverages = 1, numPeriodGrouping = 1, kargs...)
  return SimplePreProcessing(freq, frames, numAverages, loadasreal, spectralLeakageCorrection,
   bgCorrection, numPeriodAverages, numPeriodGrouping, algo, Channel{Any}(32))
end

take!(algo::SimplePreProcessing) = take!(algo.outputChannel)
function put!(algo::SimplePreProcessing, f::MPIFile)
  u = getMeasurementsFD(f, frequencies=algo.freq, frames=algo.frames, numAverages=algo.numAverages,
                          loadasreal=algo.loadasreal, spectralLeakageCorrection=algo.spectralLeakageCorrection,
                          bgCorrection=algo.bgCorrectionInternal, numPeriodAverages=algo.numPeriodAverages, numPeriodGrouping=algo.numPeriodGrouping)
  put!(algo.algo, u)
  put!(algo.outputChannel, take!(algo.algo))
end
function put!(algo::SimplePreProcessing, data::Array)
  # TODO do preprocessing as in MPIFile case
  processedData = foo(data)
  put!(algo.algo, processedData)
  put!(algo.outputChannel, take!(algo.algo))
end
outputType(algo::SimplePreProcessing) = IntermediateOutput()

mutable struct SimpleBackgroundCorrection <: AbstractMPIReconstructionAlgorithm
  uEmpty
  uEmptyPost
  bgFrames
  bgFramesPost
  # Pipeline
  algo::AbstractMPIReconstructionAlgorithm
  outputChannel::Channel{Any}  
end

take!(algo::SimpleBackgroundCorrection) = take!(algo.outputChannel)
function put!(algo::SimpleBackgroundCorrection, u)
  if bgFramesPost == nothing
    u = u .- uEmpty
  else
    for l=1:length(partframes)
      alpha = (partframes[l] - mean(bgFrames)) / (mean(bgFramesPost) - mean(bgFrames))
      u[:,:,l] .-=  (1-alpha).*uEmpty[:,:,1] .+ alpha.*uEmptyPost[:,:,1]
    end
  end
  put!(algo.algo, u)
  put!(algo.outputChannel, take!(algo.algo))
end

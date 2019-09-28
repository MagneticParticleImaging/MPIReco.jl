export reconstructionPeriodicMotion

function reconstructionPeriodicMotion(bSF::MPIFile, bMeas::MPIFile;
  minFreq=0, maxFreq=1.25e6, SNRThresh=-1,maxMixingOrder=-1, numUsedFreqs=-1, sortBySNR=false, recChannels=1:numReceivers(bMeas), kargs...)

  freq = filterFrequencies(bSF,minFreq=minFreq, maxFreq=maxFreq,recChannels=recChannels, SNRThresh=SNRThresh, numUsedFreqs=numUsedFreqs, sortBySNR=sortBySNR)

  @debug "selecting $(length(freq)) frequencies"

  return reconstructionPeriodicMotion(bSF, bMeas, freq; kargs...)
end


"""
    reconstructionPeriodicMotion(bSF::MPIFile, bMeas::MPIFile, freq::Array{Int64,1};
				bEmpty=nothing, frBG=nothing,
				alpha::Float64=3.0, choosePeak::Int64=1, frames::Int64=1,
				samplingPrecision::Bool=true, windowType::Int64=1,
				bSFFrequencyAnalysis::MPIFile=bSF,higherHarmonic::Int64=1,
				kargs...)

	Performs multi-patch reconstruction of raw data from an object with periodic motion

- bSF:			System functions for reconstruction
- bMeas:		Raw data of the measurement
- freq:                 Selected frequencies for reconstruction
- bEmpty:		Background measurement
- bgFrames:			Background frames
- choosePeak:		Number of chosen peak for motion frequency
- alpha:        Window width relative to DF cycle
- frames: 		Selected frame
- lambda:		Regularization parameter for reconstruction
- iterations:		Number of iterations
- samplingPrecision:    true: rounding motion period to sampling precision, false: rounding to DF period precision
- windowFunction:       1: Hann window, 2:FT1A05, 3: Rectangle
- bSFFrequencyAnalysis: System function for frequency analysis

"""
function reconstructionPeriodicMotion(bSF::MPIFile, bMeas::MPIFile, freq::Array{Int64,1};
				bEmpty=nothing, bgFrames=nothing,
				alpha::Float64=3.0, choosePeak::Int64=1,
                frames::UnitRange=1:acqNumFrames(bMeas),
				samplingPrecision::Bool=true, windowType::Int64=1,
				bSFFrequencyAnalysis::MPIFile=bSF, higherHarmonic::Int64=1,
				kargs...)

  FFP = squeeze(acqOffsetFieldShift(bMeas))

  motFreq = getMotionFreq(bSFFrequencyAnalysis, bMeas, choosePeak) ./ higherHarmonic
  tmot = getRepetitionsOfSameState(bMeas, motFreq, frames)

  # sort measured data in virtual frames
  uReco = getMeasurementsMotionCompFD(bMeas, motFreq, tmot, freq, frames, alpha,
                                      samplingPrecision, windowType)

  # subtract background measurement
  if bEmpty != nothing
    uEmpty = getMeasurementsFD(bEmpty, frequencies=freq, frames=1, numAverages=1, spectralLeakageCorrection=true)
    if bgFrames == nothing
      numFrames = acqNumPeriodsPerPatch(bMeas)
      bgFrames = [1+(i-1)*numFrames:i*numFrames for i=1:acqNumPatches(bBG)]
    end
    for i=1:acqNumPatches(bMeas)
      uReco[:,i,:] = uReco[:,i,:] .- mean(uEmpty[:,bgFrames[i],:], dims=2)
    end
  end

  P = numDFPeriodsInMotionCycle(motFreq, frames, dfCycle(bMeas))

  mapping = collect(1:acqNumPatches(bMeas))
  resortedInd = zeros(Int64, acqNumPatches(bMeas), P)
  for i=1:acqNumPatches(bMeas)
    resortedInd[i,:] = unflattenOffsetFieldShift(FFP)[i][1:P]
  end
  FFOp = MultiPatchOperator(bSF, bMeas, freq, false,
			    indFFPos=resortedInd[:,1],
			    FFPos=FFP[:,resortedInd[:,1]], mapping=mapping,
			    FFPosSF=FFP[:,resortedInd[:,1]])

  image = initImage(bSF[1],bMeas,size(uReco,3),acqNumAverages(bMeas),FFOp.grid,false)
  c =  reconstruction(FFOp, uReco; kargs...)
  writePartToImage!(image, c, 1, 1:size(uReco,3), acqNumAverages(bMeas))

  return image
end

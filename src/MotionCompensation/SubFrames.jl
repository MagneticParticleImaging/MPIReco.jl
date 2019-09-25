export currentFramePatchPeriod, getRepetitionsOfSameState, manualCorrection, FT1A05,
       determineWindow, getMeasurementsMotionCompFD, getPeriod

"""
    currentFramePatchPeriod(currentTimeInms, numPeriodsPerFrame, numPeriodsPerPatch,
                            numPatches, tsc[; DFCyclesSwitch==7])
Determine frame, patch, and period from time within the measurement sequence.
"""
function currentFramePatchPeriod(currentTimeInms, numPeriodsPerFrame, numPeriodsPerPatch,
                                 numPatches, tsc; DFCyclesSwitch=7)

  fr = ceil(Int, currentTimeInms/((numPeriodsPerFrame+DFCyclesSwitch*numPatches)*tsc))
  rest = currentTimeInms-(fr-1)*(numPeriodsPerFrame+DFCyclesSwitch*numPatches)*tsc
  patch = ceil(Int, rest/((numPeriodsPerPatch+DFCyclesSwitch)*tsc))
  rest2 = rest-(patch-1)*(numPeriodsPerPatch+DFCyclesSwitch)*tsc
  period = rest2/tsc
  if currentTimeInms == 0
    fr = 1
    patch = 1
    period = 1
  end

  return fr, patch, period
end

"""
    calculateRepetitionsInFramePatchPeriod(totalduration, motFreq, currentTimeInms,
                      endtime, numPeriodsPerPatch, numPeriodsPerFrame, numPatches,tsc)

Based on the starting point (currentTimeInms) in time, the repetitions are
calculated from the motion frequency at that moment (entries in motFreq) until
the endtime is reached. The times are converted to frame, patch, and period.

"""
function calculateRepetitionsInFramePatchPeriod(totalduration, motFreq, currentTimeInms, endtime,
                                        numPeriodsPerPatch, numPeriodsPerFrame, numPatches, tsc)

  tmot = zeros(ceil(Int,totalduration*maximum(motFreq))+1,3)
  j = 1
  while (currentTimeInms < endtime)
    (frame,patch,period) = currentFramePatchPeriod(currentTimeInms, numPeriodsPerFrame,
                                            numPeriodsPerPatch, numPatches,tsc)
    currentTimeInms += 1 / motFreq[frame,patch]

     if (period+1/(motFreq[frame,patch])/tsc < numPeriodsPerPatch-1 && period >= 2)
       tmot[j,1] = frame
       tmot[j,2] = patch
       tmot[j,3] = period
       j += 1
     end
  end
  return tmot
end

function calculateBeginningInms(firstFrame, tsc, numPeriodsPerPatch, numPatches, DFCyclesSwitch)
  return (firstFrame-1)*tsc*(numPeriodsPerPatch*numPatches+DFCyclesSwitch*numPatches)
end

"""
    getRepetitionsOfSameState(bMeas::MPIFile,motFreq,firstFrame,lastFrame[;DFCyclesSwitch==7])

Based on the motion frequency the repetitions of the same state are calculated. I.e. the frame, patch, and period of the repetition

#Arguments:
- `motFreq::Float`: Assumed motion frequency in Hz
"""
function getRepetitionsOfSameState(bMeas::MPIFile, motFreq, firstFrame, lastFrame;
                                   DFCyclesSwitch=7)

  numPeriodsPerFrame = acqNumPeriodsPerFrame(bMeas)
  numPatches = acqNumPatches(bMeas)
  numPeriodsPerPatch = acqNumPeriodsPerPatch(bMeas)
  endtime = (lastFrame)*dfCycle(bMeas)*(numPeriodsPerPatch*numPatches+DFCyclesSwitch*numPatches)
  currentTimeInms = calculateBeginningInms(firstFrame, dfCycle(bMeas), numPeriodsPerPatch,
                                           numPatches, DFCyclesSwitch)
  totalduration = endtime-currentTimeInms

  tmot = calculateRepetitionsInFramePatchPeriod(totalduration, motFreq, currentTimeInms,
                          endtime, numPeriodsPerPatch, numPeriodsPerFrame, numPatches,dfCycle(bMeas))

  return tmot
end


function FT1A05(windowLength)
  window = zeros(windowLength)
  n = 0:windowLength-1
  window = 0.2769-0.5261.*cos.(2*pi.*n./windowLength)+0.1971.*cos.(4*pi.*n./windowLength)
  return window
end

"""
    determineWindow(totalLength, filledLength, window)

Calculate window function with different window types of length filledLength.
If totalLength > filledLength, window function is filled with zeros.
"""
function determineWindow(totalLength, filledLength, windowType)
  window = zeros(totalLength)
  empty = 0#floor(Int,(totalLength-filledLength)/2)
  if windowType == 1
    window[empty+1:empty+filledLength] = hanning(filledLength)
  elseif windowType == 2
    window[empty+1:empty+filledLength] = FT1A05(filledLength)
  elseif windowType == 3
    window[empty+1:empty+filledLength] = ones(filledLength)
  end
  return window
end

"""
    numDFPeriodsInMotionCycle(motFreq, firstFrame, lastFrame, tsc)

Based on the motion frequency in every frame and patch, the maximum number of
complete DF cycles is caluclated for the completion of one motion cycle

"""
function numDFPeriodsInMotionCycle(motFreq, firstFrame, lastFrame, tsc)
  return floor(Int,1/minimum(motFreq[firstFrame:lastFrame,:])/tsc)
end

"""
    loadingDataIfNecessary!(bMeas::MPIFile,ui,oldFrame,currentFrame)
The data are loaded frame-wise
"""
function loadingDataIfNecessary!(bMeas::MPIFile, ui, oldFrame, currentFrame)

  if (currentFrame != oldFrame)
    ui = getMeasurements(bMeas,frames=Int(currentFrame),spectralLeakageCorrection=false)
    ui = reshape(ui,size(ui)[1],3,Int(size(ui)[3]/acqNumPatches(bMeas)),acqNumPatches(bMeas))
    oldFrame = currentFrame
  end

  return nothing#oldFrame, ui
end

"""
   calculateSubDFPeriodShift(currentPeriod, incr, samplesInOneDFCycle)

currentPeriod is given as floating number (not an Integer although this seems to be right).
The decimal digits give the shift < DF duration.
"""
function calculateSubDFPeriodShift(currentPeriod, incr, samplesInOneDFCycle, samplingPrecision)
  delta = floor(Int,(currentPeriod+incr-floor(Int,currentPeriod+incr))*samplesInOneDFCycle)
  if samplingPrecision == false
    delta = 0
  end
  return delta
end

function fillDataIntoVirtualFrame!(ui, ufinal, numMotPeriods, incrementPerPeriod,
                   samplesInOneDFCycle, currentPeriod, currentPatch, windowFunction,
                   count, samplingPrecision)

  for period = 1:numMotPeriods
    incr = incrementPerPeriod*(period-1)
    delta = calculateSubDFPeriodShift(currentPeriod,incr,samplesInOneDFCycle,samplingPrecision)

    centralPeriod = floor(Int,currentPeriod+incr)
    mixture = incr - centralPeriod

    buf = cat(ui[:,:,centralPeriod-1,currentPatch],
              ui[:,:,centralPeriod,currentPatch],
              ui[:,:,centralPeriod+1,currentPatch],
              ui[:,:,centralPeriod+2,currentPatch],dims=1)

    buf2 = broadcast(*,buf[delta+1:delta+size(ui,1)*3,:,:,:],windowFunction)
    buf3 = circshift(buf2,(delta,0,0,0))

    ufinal[:,:,period,currentPatch] += buf3[1:size(ui,1),:,:,:] +
                                       buf3[size(ui,1)+1:2*size(ui,1),:,:,:] +
                                       buf3[2*size(ui,1)+1:end,:,:,:]
  end
  return nothing
end

"""
    normalizeSignalLevelandFormatData(ufinal,numPatches,freq,count)

Normalize data based on number of repetitions to ensure same signal level for all patches and reformatting
"""
function normalizeSignalLevelandFormatData(ufinal, numPatches, freq, count)
  for patch in 1:numPatches
    ufinal[:,:,:,patch] /= count[patch]
  end
  ufinal = rfft(ufinal,1)
  ufinal = permutedims(reshape(ufinal,size(ufinal)[1]*size(ufinal)[2],
                               size(ufinal)[3],size(ufinal)[4])[freq,:,:],[1,3,2])
  return ufinal
end

"""
	getMeasurementsMotionCompFD(bMeas::MPIFile, motFreq, tmot, freq, firstFrame, lastFrame, sigma,
                      samplingPrecision, windowType)
	Averages over raw data corresponding to the same motion state and creates the virtual frames

- bMeas:		Meas data
- motFreq:		Array containing dominant motion frequencies for each patch and frame
- tmot:			Array containing repetitions of the same state (frame,patch,period)
- freq:			Selected frequencies for reconstruction
- firstFrame/lastframe:	Selected frames
- Δt:			Window width for spectral leakage correction (Δt = 3 <=> window width = 3*DF repetition time)
- samplingPrecision:	true: rounding motion period to sampling precision, false: rounding to DF period precision
- window:		1: Hann window, 2:FT1A05, 3: Rectangle

"""
function getMeasurementsMotionCompFD(bMeas::MPIFile, motFreq, tmot, freq, firstFrame, lastFrame, Δt,
                           samplingPrecision, windowType)
  oldFrame = firstFrame

  numMotPeriods = numDFPeriodsInMotionCycle(motFreq, firstFrame, lastFrame, dfCycle(bMeas))

  ui = getMeasurements(bMeas,frames=oldFrame,spectralLeakageCorrection=false)
  ui = reshape(ui,size(ui)[1],3,Int(size(ui)[3]/acqNumPatches(bMeas)),acqNumPatches(bMeas))
  ufinal = zeros(Float32,size(ui)[1],size(ui)[2],numMotPeriods,size(ui)[4])
  count = zeros(acqNumPatches(bMeas))
  window = determineWindow(size(ui)[1]*3,floor(Int,Δt*size(ui)[1]), windowType)

  #numPeriods = acqNumPeriodsPerPatch(b)

  # collect data from all repetitions from first till last frame and fill into virtual frame ufinal
  for i = 1:size(tmot)[1]
    currentFrame=Int(tmot[i,1])
    if currentFrame > firstFrame-1 && currentFrame < lastFrame+1
      loadingDataIfNecessary!(bMeas,ui,oldFrame,currentFrame)
      #(oldFrame, ui) = loadingDataIfNecessary(bMeas,oldFrame,currentFrame,ui)

      currentPatch = Int(tmot[i,2])
      # Counting is required for averaging to ensure same signal level for all patches
      count[currentPatch] += 1
      currentMotFreq = motFreq[currentFrame,currentPatch]
      # The number of full DF cycles within one motion cycle varies for different patches and frames
      # For consistent images, the same number is required for all patches. Thus for higher motion in a patch, the increment per Period is lowered
      incrementPerPeriod = 1/numMotPeriods*(1/currentMotFreq/dfCycle(bMeas))
      currentPeriod = tmot[i,3]
      fillDataIntoVirtualFrame!(ui, ufinal, numMotPeriods, incrementPerPeriod,
                               size(ui)[1], currentPeriod, currentPatch, window,
                               count, samplingPrecision)
    end
  end

  return normalizeSignalLevelandFormatData(ufinal, acqNumPatches(bMeas), freq, count)
end

function getPeriod(bMeas::MPIFile, freq)
  period = 1 / (dfCycle(bMeas)*freq)
  return period
end

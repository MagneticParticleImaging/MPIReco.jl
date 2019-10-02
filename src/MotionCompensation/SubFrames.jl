export currentFramePatchPeriod, getRepetitionsOfSameState, manualCorrection, FT1A05,
       determineWindow, getMeasurementsMotionCompFD, getPeriod


struct LocationInMeasurement
  frame::Int
  patch::Int
  period::Float64
end

"""
    currentFramePatchPeriod(currentTimeInms, numPeriodsPerFrame, numPeriodsPerPatch,
                            numPatches, tsc[; DFCyclesSwitch==7])
Determine frame, patch, and period from time within the measurement sequence.
"""
function currentFramePatchPeriod(currentTimeInms, numPeriodsPerFrame, numPeriodsPerPatch,
                                 numPatches, tsc; DFCyclesSwitch=7)

  frame = ceil(Int, currentTimeInms/((numPeriodsPerFrame+DFCyclesSwitch*numPatches)*tsc))
  rest = currentTimeInms-(frame-1)*(numPeriodsPerFrame+DFCyclesSwitch*numPatches)*tsc
  patch = ceil(Int, rest/((numPeriodsPerPatch+DFCyclesSwitch)*tsc))
  rest2 = rest-(patch-1)*(numPeriodsPerPatch+DFCyclesSwitch)*tsc
  period = rest2/tsc
  if currentTimeInms == 0
    frame = 1
    patch = 1
    period = 1.0
  end
  return LocationInMeasurement(frame, patch, period)
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

  #N = ceil(Int, totalduration*maximum(motFreq)) + 1
  tmot = Vector{LocationInMeasurement}()

  while (currentTimeInms < endtime)
    loc = currentFramePatchPeriod(currentTimeInms, numPeriodsPerFrame,
                                            numPeriodsPerPatch, numPatches,tsc)

    currentTimeInms += 1 / motFreq[loc.frame,loc.patch]

     if (loc.period+1/(motFreq[loc.frame,loc.patch])/tsc < numPeriodsPerPatch-1 &&
         loc.period >= 2)
       push!(tmot, loc)
     end
  end
  return tmot
end

function calculateBeginningInms(firstFrame, tsc, numPeriodsPerPatch, numPatches, DFCyclesSwitch)
  return (firstFrame-1)*tsc*(numPeriodsPerPatch*numPatches+DFCyclesSwitch*numPatches)
end

"""
    getRepetitionsOfSameState(bMeas::MPIFile,motFreq,frames[;DFCyclesSwitch==7])

Based on the motion frequency the repetitions of the same state are calculated. I.e. the frame, patch, and period of the repetition

#Arguments:
- `motFreq::Float`: Assumed motion frequency in Hz
"""
function getRepetitionsOfSameState(bMeas::MPIFile, motFreq, frames::UnitRange;
                                   DFCyclesSwitch=7)

  numPeriodsPerFrame = acqNumPeriodsPerFrame(bMeas)
  numPatches = acqNumPatches(bMeas)
  numPeriodsPerPatch = acqNumPeriodsPerPatch(bMeas)
  endtime = last(frames)*dfCycle(bMeas)*(numPeriodsPerPatch*numPatches+DFCyclesSwitch*numPatches)
  currentTimeInms = calculateBeginningInms(first(frames), dfCycle(bMeas), numPeriodsPerPatch,
                                           numPatches, DFCyclesSwitch)
  totalduration = endtime-currentTimeInms

  tmot = calculateRepetitionsInFramePatchPeriod(totalduration, motFreq, currentTimeInms,
                          endtime, numPeriodsPerPatch, numPeriodsPerFrame, numPatches, dfCycle(bMeas))

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
    numDFPeriodsInMotionCycle(motFreq, frames, tsc)

Based on the motion frequency in every frame and patch, the maximum number of
complete DF cycles is caluclated for the completion of one motion cycle

"""
function numDFPeriodsInMotionCycle(motFreq::Matrix, frames, tsc)
  return floor(Int, 1.0 / (minimum(motFreq[frames,:])*tsc) )
end

function numDFPeriodsInMotionCycle(motFreq::Number, tsc)
  return floor(Int, 1.0 / (motFreq*tsc) )
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

function fillDataIntoVirtualFrame!(ui, ufinal, normalization, numMotPeriods,
                   incrementPerPeriod, samplesInOneDFCycle, currentPeriod,
                   currentPatch, window, samplingPrecision)

  numDFCyclesPerWindowLength = div(length(window), size(ui,1))
  N = size(ui,1)

  for period = 1:numMotPeriods
    incr = incrementPerPeriod*(period-1)
    delta = calculateSubDFPeriodShift(currentPeriod,incr,samplesInOneDFCycle,samplingPrecision)

    centralPeriod = floor(Int, currentPeriod+incr)
    shiftedWindow = circshift(window, delta)

    for recChan = 1:size(ui,2)
      p = (centralPeriod-1):(centralPeriod-1+numDFCyclesPerWindowLength)
      buf = vec(ui[:,recChan,p,currentPatch])
      buf2 = circshift(buf[delta+1:delta+N*numDFCyclesPerWindowLength,:], delta)
      buf3 = buf2 .* shiftedWindow

      for l=1:numDFCyclesPerWindowLength
        ufinal[:,recChan,period,currentPatch] += buf3[((l-1)*N+1):(l*N)]
        normalization[:,recChan,period,currentPatch] +=
                             shiftedWindow[((l-1)*N+1):(l*N)]
      end
    end
  end
  return nothing
end

"""
    normalizeSignalLevelandFormatData(ufinal,numPatches,freq,count)

Normalize data based on number of repetitions to ensure same signal level for all patches and reformatting
"""
function normalizeSignalLevelandFormatData(ufinal, normalization, freq)
  if any(normalization .== 0)
    error("Could not construct virtual frame!")
  end 
  ufinal ./= normalization

  ufinal = rfft(ufinal,1)
  ufinal = permutedims(reshape(ufinal,size(ufinal)[1]*size(ufinal)[2],
                               size(ufinal)[3],size(ufinal)[4])[freq,:,:],[1,3,2])
  return ufinal
end

"""
	getMeasurementsMotionCompFD(bMeas::MPIFile, motFreq, tmot, freq, frames, sigma,
                      samplingPrecision, windowType)
	Averages over raw data corresponding to the same motion state and creates the virtual frames

- bMeas:		Meas data
- motFreq:		Array containing dominant motion frequencies for each patch and frame
- tmot:			Array containing repetitions of the same state (frame,patch,period)
- freq:			Selected frequencies for reconstruction
- frames:       Selected frames
- alpha:        Window width relative to the DF cycle.
- samplingPrecision:	true: rounding motion period to sampling precision, false: rounding to DF period precision
- window:		1: Hann window, 2:FT1A05, 3: Rectangle

"""
function getMeasurementsMotionCompFD(bMeas::MPIFile, motFreq, tmot, freq, frames::UnitRange,
                                     alpha, samplingPrecision, windowType)
  oldFrame = first(frames)

  numMotPeriods = numDFPeriodsInMotionCycle(motFreq, frames, dfCycle(bMeas))

  ui_ = getMeasurements(bMeas, frames=first(frames), spectralLeakageCorrection=false)
  ui = reshape(ui_, size(ui_,1), 3, div(size(ui_,3), acqNumPatches(bMeas)), acqNumPatches(bMeas))
  ufinal = zeros(Float32, size(ui,1), size(ui,2), numMotPeriods, size(ui,4))
  normalization = zeros(Float32, size(ui,1), size(ui,2), numMotPeriods, size(ui,4))
  count = zeros(acqNumPatches(bMeas))
  numDFCyclesPerWindowLength = ceil(Int, alpha)
  window = determineWindow(size(ui,1)*numDFCyclesPerWindowLength,
                           floor(Int, alpha*size(ui,1)), windowType)

  # collect data from all repetitions from first till last frame and fill into virtual frame ufinal
  for i = 1:length(tmot)
    currentFrame = tmot[i].frame
    if currentFrame > first(frames)-1 && currentFrame < last(frames)+1
      # loading new data if frame has changed
      if currentFrame != oldFrame
        ui[:] = vec(getMeasurements(bMeas, frames=currentFrame, spectralLeakageCorrection=false))
        oldFrame = currentFrame
      end

      currentPatch = tmot[i].patch

      currentMotFreq = motFreq[currentFrame,currentPatch]
      # The number of full DF cycles within one motion cycle varies for different patches and frames
      # For consistent images, the same number is required for all patches.
      # Thus for higher motion in a patch, the increment per Period is lowered
      incrementPerPeriod = numDFPeriodsInMotionCycle(currentMotFreq, dfCycle(bMeas)) / numMotPeriods

      currentPeriod = tmot[i].period
      fillDataIntoVirtualFrame!(ui, ufinal, normalization, numMotPeriods, incrementPerPeriod,
                               size(ui,1), currentPeriod, currentPatch, window,
                               samplingPrecision)
    end
  end

  return normalizeSignalLevelandFormatData(ufinal, normalization, freq)
end

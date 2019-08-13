export reconstructionPeriodicMotion

"""
	reconstructionPeriodicMotion(b::MPIFile, bSF::MPIFile, bBG::MPIFile, frBG::Array{UnitRange{Int64},1}, choosePeak::Int64, sigma::Float64, freq::Array{Int64,1},recoFrame::Int64;lambd=0.1,iterations=2,samplingPrecision=true,windowFunction=1,bSFFrequencyAnalysis=bSF)

	Performs multi-patch reconstruction of raw data from an object with periodic motion

- bMeas:		Raw data of the measurement
- bSF:			System functions for reconstruction
- bBG:			Background measurement
- frBG:			Background frames
- choosePeak:		Number of chosen peak for motion frequency
- alpha:                Window width for spectral leakage correction alpha = 1 <=> 3*DF repetition time
- freq:                 Selected frequencies for reconstruction
- recoFrame: 		Selected frame
- lambda:		Regularization parameter for reconstruction
- iterations:		Number of iterations
- samplingPrecision:    true: rounding motion period to sampling precision, false: rounding to DF period precision
- windowFunction:       1: Hann window, 2:FT1A05, 3: Rectangle+
- bSFFrequencyAnalysis: System function for frequency analysis

"""
function reconstructionPeriodicMotion(bMeas::MPIFile, bSF::MPIFile, bBG::MPIFile,
	                    	frBG::Array{UnitRange{Int64},1}, choosePeak::Int64,
				alpha::Float64, freq::Array{Int64,1},recoFrame::Int64;
				lambda=0.1, iterations=2, 
				samplingPrecision=true, windowType=1,
				bSFFrequencyAnalysis=bSF,higherHarmonic=1)

  FFP = squeeze(acqOffsetFieldShift(bMeas))

  uBG = getMeasurementsFD(bBG, frequencies=freq, frames=1, numAverages=1, spectralLeakageCorrection=true)

  motFreq = getMotionFreq(bMeas,bSFFrequencyAnalysis,choosePeak)./higherHarmonic
  tmot = getRepetitionsOfSameState(bMeas,motFreq,recoFrame,recoFrame)

  # sort measured data in virtual frames
  resortedInd = zeros(Int64,acqNumPatches(bMeas),floor(Int,getPeriod(bMeas, motFreq[1,1])))
  for i=1:acqNumPatches(bMeas)
    resortedInd[i,:] = unflattenOffsetFieldShift(FFP)[i][1:floor(Int,getPeriod(bMeas, motFreq[1,1]))]
  end
  uReco = getMeasurementsMotionCompFD(bMeas, motFreq, tmot, freq, recoFrame, recoFrame, alpha,
                            samplingPrecision, windowType)
  #println(size(uReco))
  for i=1:acqNumPatches(bMeas)
    uReco[:,i,:] = uReco[:,i,:] .- mean(uBG[:,frBG[i],:], dims=2)
  end
  mapping = collect(1:acqNumPatches(bMeas)) #ones(Int,acqNumPatches(bMeas))
  FFOp = MultiPatchOperator(bSF, bMeas, freq, false,
			    indFFPos=resortedInd[:,1],
			    FFPos=FFP[:,resortedInd[:,1]], mapping=mapping, 
			    FFPosSF=FFP[:,resortedInd[:,1]])

  c_ =  reconstruction(FFOp, uReco, λ=lambda, iterations=iterations)
  return reshape(c_, shape(FFOp.grid)..., :)
end

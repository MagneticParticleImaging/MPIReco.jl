export reconstructionPeriodicMotion

"""
	reconstructionPeriodicMotion(b::MPIFile, bSF::MPIFile, bBG::MPIFile, frBG::Array{UnitRange{Int64},1}, choosePeak::Int64, sigma::Float64, freq::Array{Int64,1},recoFrame::Int64;lambd=0.1,iterations=2,samplingPrecision=true,windowFunction=1,bSFFrequencyAnalysis=bSF)

	Performs multi-patch reconstruction of raw data from an object with periodic motion

- b:			Raw data
- bSF:			System function for reconstruction
- bBG:			Background measurement
- frBG:			Background frames
- choosePeak:		Number of chosen peak for motion frequency
- sigma:                Window width for spectral leakage correction sigma = 1 <=> 3*DF repetition time
- freq:                 Selected frequencies for reconstruction
- recoFrame: 		Selected frame
- lambd:		Regularization parameter for reconstruction
- iterations:		Number of iterations
- samplingPrecision:    true: rounding motion period to sampling precision, false: rounding to DF period precision
- windowFunction:       1: Hann window, 2:FT1A05, 3: Rectangle+
- bSFFrequencyAnalysis: System function for frequency analysis

"""
function reconstructionPeriodicMotion(b::MPIFile, bSF::MPIFile, bBG::MPIFile,
	                    frBG::Array{UnitRange{Int64},1}, choosePeak::Int64,
						sigma::Float64, freq::Array{Int64,1},recoFrame::Int64;
						lambda=0.1, iterations=2, samplingPrecision=true, windowType=1,
						bSFFrequencyAnalysis=bSF)

  FFP = ffPos(b)

  uBG = getMeasurementsFD(bBG, frequencies=freq, frames=1, numAverages=1, spectralLeakageCorrection=true)

  motFreq = getMotionFreq(b,bSFFrequencyAnalysis,choosePeak)
  tmot = getRepetitionsOfSameState(motFreq,b,recoFrame,recoFrame)

  resortedInd = unflattenOffsetFieldShift(FFP)[:,1:floor(Int,getPeriod(b, motFreq[1,1]))]
  uReco = getavrgusubPeriod(motFreq, tmot, b, freq, recoFrame, recoFrame, sigma,
                            samplingPrecision, windowType)
  #println(size(uReco))
  for i=1:acqNumPatches(b)
    uReco[:,i,:] = uReco[:,i,:] .- mean(uBG[:,frBG[i],:], dims=2)
  end
  mapping = ones(Int,acqNumPatches(b))
  FFOp = FFOperator(bSF, b, freq, false,
                    #OverscanSF=[0.0179,0.01,0.0105]./2,   # not used anymore
                    #OffsetFF=[0.0179,0.01,0.0105]./2,     # not used anymore
					indFFPos=resortedInd[:,1],
					FFPos=FFP[:,resortedInd[:,1]], mapping=mapping, FFPosSF=FFP[:,resortedInd[:,1]])

 c_ =  reconstruction(FFOp, uReco, λ=lambda, iterations=iterations)
 return reshape(c_, shape(FFOp.grid)..., :)
end

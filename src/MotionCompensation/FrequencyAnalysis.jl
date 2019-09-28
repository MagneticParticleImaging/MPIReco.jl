export findComponentsWithSignifFreq, getMotionFreq, findPeaks,
       bins2hertz, hertz2bins, generateSpectrum,
       parabelFit, parabelMax

function findComponentsWithSignifFreq(u)
  uFD = zeros(size(u,1),Int64(size(u,2)/2))
  for i = 1:size(u,1)
    uFD[i,:] = generateSpectrum(squeeze(u[i,:]))
  end
  value = zeros(size(uFD,1))
  for i = 1:size(uFD,1)
    value[i] = maximum(squeeze(uFD[i,:]))
  end
  value[1] = 0
  return findmax(value)
end

"""
    getMotionFreq(bSF::MPIFile, bMeas::MPIFile, choosePeak::Int)

Determines frequency [Hz] of signal modulation for each patch and frame

...
- `bSF::MPIFile`: Systemfunction ==> required because the motion frequency is only determined from components with high SNR
- `bMeas::MPIFile`: measurement
- `choosePeak::Int`: Within the measurement signal different modulations can be found and also higher harmonics.
 Choose the number of the peak you want
...
"""
function getMotionFreq(bSF::MPIFile, bMeas::MPIFile, choosePeak::Int,
		       SNRThreshMF = 10, minFreqMF = 80e3, recChannelsMF = [2]) # frequency selection

    freqs = filterFrequencies(bSF,SNRThresh=SNRThreshMF,minFreq=minFreqMF,recChannels=recChannelsMF)
    repFreq = 1 / dfCycle(bMeas);
    peaksInHz = 0
    selectedPeakInHz = zeros(acqNumFrames(bMeas),acqNumPatches(bMeas))
    uspec = getMeasurementsFD(bMeas, frequencies=freqs, numAverages=1, spectralLeakageCorrection=false)
    norm, ind = findComponentsWithSignifFreq(squeeze(uspec[:,1:acqNumPeriodsPerPatch(bMeas),1]))
    @debug "Maximum freq" ind

	for j=1:acqNumFrames(bMeas)
    	for i = 1:acqNumPatches(bMeas)
        	uFT = generateSpectrum(squeeze(uspec[ind,(i-1)*acqNumPeriodsPerPatch(bMeas)+1:i*acqNumPeriodsPerPatch(bMeas),j]))
        	frind = collect(1:acqNumPeriodsPerPatch(bMeas)/2)
       		frHz = bins2hertz(frind,repFreq,acqNumPeriodsPerPatch(bMeas))

            peaksInInt = findPeaks(uFT,stdTh=1)
            peaksInHz = bins2hertz(peaksInInt, 1/dfCycle(bMeas), acqNumPeriodsPerPatch(bMeas))
            @debug "Found peaks: " peaksInHz

			selectedPeak = parabelMax(parabelFit(peaksInInt[choosePeak]-1:peaksInInt[choosePeak]+1,log.(uFT[peaksInInt[choosePeak]-1:peaksInInt[choosePeak]+1])))
			selectedPeakInHz[j,i] = bins2hertz(selectedPeak,repFreq,acqNumPeriodsPerPatch(bMeas))
		end
    end
    return selectedPeakInHz
end

function parabelFit(x,y)
    p = [x.^2 x ones(x)]\y
    return p
end

function parabelMax(par)
    return -par[2]/2/par[1]
end

function generateSpectrum(data)
  data = applyHannWindow(broadcast(-, data, mean(data)))
  ftc = fft(data)

  ft = ftc
  ft = abs.(ftc[2:div(end,2)+1])
  helpv = size(ft,1)
  ft[1:end] += abs.(ftc[end:-1:end-helpv+1])

  return ft
end

function applyHannWindow(data)

    N = size(data,1)
    n = (0:N-1)

    return data .* (0.5*( 1 .-cos.(2*pi*n/(N-1)) ))
end


function reduceHighSpectralValuesToPeaks(data, isHigh)
  diracpeaks = []
  idx = 1
  while !isempty(findall(isHigh[idx:end]))
      from = findfirst(isHigh[idx:end]) + idx - 1
      till = findfirst(isHigh[from+1:end].==0) + from - 1
      push!(diracpeaks, findmax(data[from:till])[2]+from-1)
      idx = till + 1
  end
  return diracpeaks
end

function enforceMinDistanceOfPeaks(diracpeaks, minBinDist, data)
  keep = false
  while sum(.!keep) > 0
      keep = ( diracpeaks[2:end]-diracpeaks[1:end-1] ) .> minBinDist
      push!(keep,true) # append last index
      # keep the greater ones
      rem = findall(.!keep)
      @debug rem
      for i in rem
	  @debug diracpeaks[i]
          if data[diracpeaks[i]] > data[diracpeaks[i+1]]
              keep[i] = true
              keep[i+1] = false
          end
      end
      diracpeaks = diracpeaks[keep .== true]

  end
  return diracpeaks
end

"""
	findPeaks(data; stdTh=1.5, minBinDist=0, minBin=1)
Find peaks in the data.

- data:         column vector of data
- stdTh:        peaks are at least stdTh*standarddeviation above the mean
- minBinDist:   min number of bins between two peaks
- minBin:       no peaks in bins below minBin
- return:       peak bins
"""
function findPeaks(data; stdTh=1.5, minBinDist=0, minBin=1)
    # smooth data??

    if length(findall(data[1:end-1].<data[2:end])) == 0 || length(findall(data[1:end-1].>data[2:end])) == 0
	warn("This should not happen... No variation of data")
	return []
    end

    # exclude beginning and end of data
    frst = maximum( [findfirst( data[1:end-1].<data[2:end] ), ceil(Int32, minBin)] )
    lst = minimum( [findlast( data[1:end-1].>data[2:end] ), floor(Int32, length(data)-minBin+1)] ) # last entry of data never used
    data = data[frst:lst]

    # detrending
    #x = 1:length(data)
    #linPar = [x ones(x)]\data
    #data = data - (linPar[1]*x+linPar[2])

    # mean,std
    mu = mean(data)
    sigma = std(data)

    # threshold on diff to mean
    peaks = data .> (mu + stdTh*sigma)
    if isempty(peaks)
        @error "No peaks were found!"
    end

    # no maximum at last position
    peaks[end] = false

    diracpeaks = reduceHighSpectralValuesToPeaks(data,peaks)
    diracpeaks = enforceMinDistanceOfPeaks(diracpeaks,minBinDist,data)

    diracpeaks = diracpeaks .+ frst.-1
    return diracpeaks
end


"""
	bins2hertz(bin, fs, nSamples)
Convert FFT bin index to Hertz.

- bin:      bin index
- fs:       sampling frequency
- nSamples: number of samples FFT was applied to
- return:   Hertz
"""
function bins2hertz(bin, fs, nSamples)
  return fs/nSamples*bin
end


"""
	hertz2bins(hz, fs, nSamples)
Convert Hertz to FFT bin index.

- hz:       Hertz
- fs:       sampling frequency
- nSamples: number of samples FFT was applied to
- return:   bins
"""
function hertz2bins(hz, fs, nSamples)
  return nSamples/fs*hz
end

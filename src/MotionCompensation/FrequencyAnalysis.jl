export findComponentsWithSignifFreq, getMotionFreq, findPeaks, peakHist,
       bins2hertz, hertz2bins, generateSpectrum, mpiFFT, dynComps, hannWindow,
       parabelFit, parabelMax, peakComps, findNFreqwithhighestSNRinfreqs,
       mpiFFTpixel, mpiFFTeval

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

"""
    findNFreqwithhighestSNRinfreqs(bSF::MPIFile, freqs, N)

Return indices of N frequency components with highest SNR.
"""
function findNFreqwithhighestSNRinfreqs(bSF::MPIFile, freqs, N)

  snr_freq = reshape(getSNR(bSF),size(getSNR(bSF),1)*3,1)
  indinfreqs = Int[]
  buf = 0
  for ze =1:N
    if ze == 1
      buf = maximum(snr_freq)
    end
    maxsnrcomps=indmax(snr_freq)
    for i = 1:length(freqs)
      if (freqs[i]==maxsnrcomps)
        push!(indinfreqs,i)
      end
    end
    if (ze ==N)
      @debug "" snr_freq[maxsnrcomps]/buf
    end
    snr_freq[maxsnrcomps]=0
  end

  return indinfreqs
end

function parabelFit(x,y)
    p = [x.^2 x ones(x)]\y
    return p
end

function parabelMax(par)
    return -par[2]/2/par[1]
end

"""
    dynComps(u)
Return indices of components with a variance higher than 0.1 of maximum variance
"""
function dynComps(u)
	vc = var(u,2)# ./ abs(mean(u,2));
	idx = vc .> maximum(vc)*0.1
	return find(idx)
end

"""
    peakComps(u, peak, t=0.5)
"""
function peakComps(u, peak, t=0.5)
    peak = round(Int32, peak)
	ft = mpiFFT(u, full=true, raw=true)
    c = sum(ft[:,peak:peak],2) ./ mean(ft,2)
    @debug "" size(ft[:,peak:peak])
    idx = c .> t*maximum(c)
	return find(idx)
end

function peakComps2(u, peak)
	ft = mpiFFT(u, full=true, raw=true)
    c = sum(ft[:,peak-3:peak+3]) ./ mean(ft,2)
    idx = falses(size(ft,1),1)
    for i = 1:10
        m = indmax(c)
        idx[m] = true
        c[m] = 0
    end
  return idx
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

"""
	mpiFFT(u)

Compute spectrum of temporal series of measurement vectors.

- u:        series of measurement vectors
- dc:       return DC component
- hann:     apply hann window before FFT
- return:   spectrum
"""
function mpiFFT(u; dc=false, hann=true, full=false, raw=false, comps=nothing, normalizecomps=false)
	#u = reshape(u,size(u,1)*size(u,2),size(u,3))
    if full
      idx = 1:size(u,1)
    elseif comps != nothing
      idx = comps
    else
      idx = dynComps(u)
    end
	if hann
		u = hannWindow(broadcast(-,u,mean(u,2)))
		#u = hannWindow(u-repmat(mean(u,2),1,size(u,2)))
	end
	ftc = fft(u[idx,:],2)

    ft = ftc
    ft = abs.(ftc[:,2-dc:div(end,2)+1])
    helpv = size(ft,2)-dc
    ft[:,1+dc:end] += abs.(ftc[:,end:-1:end-helpv+1])

    if (normalizecomps == true)
      maxima = maximum(ft,2)
      ft = broadcast(/,ft,maxima)
    end

    if raw
        return ft
    else
        return (mean(ft,1))'
    end
end


function mpiFFTeval(u, freqRange; dc=false, hann=true, window=63, freqdecrfac=1)

  sampleFreq = 1/(21.54e-3)
  minfreq = maximum([1, floor(Int64,hertz2bins(freqRange[1], sampleFreq, window))])
  maxfreq = ceil(Int64,hertz2bins(freqRange[2], sampleFreq, window))
  halfwindow = div(window,2)
  result = Float64[];
  freqint = floor(Int64,(maxfreq-minfreq+1)*freqdecrfac)

  f = 1
  while f < length(u)-window
    ft = mpiFFTpixel(u[f:f+window])
    m = findmax(ft[minfreq:maxfreq]) + minfreq - 1
    v = parabelMax(parabelFit(m-1:m+1,ft[m-1:m+1]))
    vhz = bins2hertz(v, sampleFreq, window)
    push!(result,vhz)
    minfreq = maximum([1,m - freqint])
    maxfreq = m + freqint
    f+=1
  end

  return result
end

function mpiFFTpixel(u; dc=false, hann=true)

  if hann
    u = hannWindow(u-mean(u))
  end

  ftc = abs.(rfft(u[:],1))

  return ftc
end

function applyHannWindow(data)

    N = size(data,1)
    n = (0:N-1)

    return data .* (0.5*( 1 .-cos.(2*pi*n/(N-1)) ))
end


"""
	hannWindow(data)
Apply a Hann window function to a vector of data. If a matrix is passed
the window is applied row-wise.

- data:	   data
- return:  windowed data as row vector(s)
"""
function hannWindow(data)
    if size(data,2) == 1
      #data = data';
      data = permutedims(data,[2,1])
    end

    N = size(data,2)
    n = (0:N-1)'

    return data .* (0.5*( 1-cos.(2*pi*n/(N-1)) ))
end


"""
	peakHist(data; stdTh=1.5,minBinDist=0,minBin=1)
Histogram of the peaks in the FFTs of each single frequency component.

- data:         freq vectors, nFreqs x nFrames
- stdTh:        peaks are at least stdTh*standarddeviation above the mean
- minBinDist:   min number of bins between two peaks
- minBin:       no peaks in bins below minBin
- return:       histogram of peaks
"""
function peakHist(data; stdTh=1.5, minBinDist=0, minBin=1)
    th = stdTh
    N = size(data,1)
    h = zeros(size(data,2))
    progmeter = Progress(N,3)
    for i = 1:N
        p = findPeaks(data[i,:]',stdTh=th,minBinDist=minBinDist,minBin=minBin)
        h[p] = h[p] + 1
	next!(progmeter)
    end
    return h
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

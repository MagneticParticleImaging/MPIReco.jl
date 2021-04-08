using MPIReco

@testset "sparse single- and multi-channel in-memory reconstruction" begin
  bSF = MPIFile("$datadir/mdf/systemMatrix.mdf")
  b = MPIFile("$datadir/mdf/measurement.mdf")
  redFactor = 0.01
  names = (:color, :x, :y, :z, :time)
  values = (1:1,
	    -21.5u"mm":1.0u"mm":21.5u"mm",
	    -21.5u"mm":1.0u"mm":21.5u"mm",
	    0.0u"mm":1.0u"mm":0.0u"mm",
	    0.0u"ms":65.28u"ms":0.0u"ms")
  valuesDF = (1:1,
	      -21.51304347826087u"mm":0.9739130434782609u"mm":-0.08695652173913043u"mm",
	      -21.51304347826087u"mm":0.9739130434782609u"mm":-0.08695652173913043u"mm",
	      -0.5u"mm":1.0u"mm":-0.5u"mm",
	      0.0u"ms":65.28u"ms":0.0u"ms")

  c1 = reconstruction(bSF, b; lambd=0.1,
		      SNRThresh = 2, frames=1:100, minFreq=80e3, nAverages=100,
		      recChannels =1:2, iterations = 3,
		      spectralLeakageCorrection=false)
  @test axisnames(c1) == names
  @test axisvalues(c1) == values
  exportImage("$imgdir/Sparse1.png", arraydata(c1[1,:,:,1,1]))

  c2 = reconstruction(bSF, b; lambd=0.1,SNRThresh = 2, frames=1:100,
		      minFreq=80e3, nAverages=100, sparseTrafo="FFT",
		      recChannels =1:2, iterations = 3,
		      spectralLeakageCorrection=false, redFactor = redFactor,
		      useDFFoV = false)
  @test axisnames(c2) == names
  @test axisvalues(c2) == values
  exportImage("$imgdir/Sparse2.png", arraydata(c2[1,:,:,1,1]))

  c3 = reconstruction(bSF, b; lambd=0.1, SNRThresh = 2, frames=1:100,
		      minFreq=80e3, nAverages=100, sparseTrafo="FFT",
		      recChannels =1:2, iterations = 3,
		      spectralLeakageCorrection=false, redFactor = redFactor,
		      useDFFoV = true)
  @test axisnames(c3) == names
  @test axisvalues(c3) == valuesDF
  exportImage("$imgdir/Sparse3.png", arraydata(c3[1,:,:,1,1]))

  c4 = reconstruction(bSF, b; lambd=0.01, SNRThresh = 3, frames=1:100,
		      minFreq=80e3, nAverages=100, sparseTrafo="DCT-IV",
		      recChannels =1:2, iterations = 1,
		      spectralLeakageCorrection=false, redFactor = redFactor,
		      useDFFoV = false)
  @test axisnames(c4) == names
  @test axisvalues(c4) == values
  exportImage("$imgdir/Sparse4.png", arraydata(c4[1,:,:,1,1]))

  c5 = reconstruction(bSF, b; lambd=0.01, SNRThresh = 3, frames=1:100,
		      minFreq=80e3, nAverages=100, sparseTrafo="DCT-IV",
		      recChannels =1:2, iterations = 1,
		      spectralLeakageCorrection=true, redFactor = redFactor,
		      useDFFoV = true)
  @test axisnames(c5) == names
  @test axisvalues(c5) == valuesDF
  exportImage("$imgdir/Sparse5.png", arraydata(c5[1,:,:,1,1]))

  c6 = reconstruction(bSF, b; lambd=0.01, SNRThresh = 3, frames=1:100,
		      minFreq=80e3, nAverages=100, sparseTrafo="DST",
		      recChannels =1:2, iterations = 1,
		      spectralLeakageCorrection=false, redFactor = redFactor,
		      useDFFoV = false)
  @test axisnames(c6) == names
  @test axisvalues(c6) == values
  exportImage("$imgdir/Sparse6.png", arraydata(c6[1,:,:,1,1]))

  c7 = reconstruction(bSF, b; lambd=0.01, SNRThresh = 3, frames=1:100,
		      minFreq=80e3, nAverages=100, sparseTrafo="DST",
		      recChannels =1:2, iterations = 1,
		      spectralLeakageCorrection=true, redFactor = redFactor,
		      useDFFoV = true)
  @test axisnames(c7) == names
  @test axisvalues(c7) == valuesDF
  exportImage("$imgdir/Sparse7.png", arraydata(c7[1,:,:,1,1]))
end

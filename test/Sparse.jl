using MPIReco
using Test
using Winston

@testset "sparse single- and multi-channel in-memory reconstruction" begin
  bSF = MPIFile("systemMatrix.mdf")
  b = MPIFile("measurement.mdf")
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
  Winston.colormap("Grays")

  c1 = reconstruction(bSF, b; lambd=0.1,
		      SNRThresh = 2, frames=1:100, minFreq=80e3, nAverages=100,
		      recChannels =1:2, iterations = 3, 
		      spectralLeakageCorrection=false)
  @test axisnames(c1) == names
  @test axisvalues(c1) == values
  p = imagesc(data(data(c1[1,:,:,1,1])), (minimum(c1),maximum(c1)))
  #display(p)

  c2 = reconstruction(bSF, b; lambd=0.1,SNRThresh = 2, frames=1:100, 
		      minFreq=80e3, nAverages=100, sparseTrafo="FFT",
		      recChannels =1:2, iterations = 3, 
		      spectralLeakageCorrection=false, redFactor = redFactor,
		      useDFFoV = false)
  @test axisnames(c2) == names
  @test axisvalues(c2) == values
  p = imagesc(data(data(c2[1,:,:,1,1])), (minimum(c2),maximum(c2)))
  #display(p)

  c3 = reconstruction(bSF, b; lambd=0.1, SNRThresh = 2, frames=1:100, 
		      minFreq=80e3, nAverages=100, sparseTrafo="FFT", 
		      recChannels =1:2, iterations = 3, 
		      spectralLeakageCorrection=false, redFactor = redFactor, 
		      useDFFoV = true)
  @test axisnames(c3) == names
  @test axisvalues(c3) == valuesDF
  p = imagesc(data(data(c3[1,:,:,1,1])), (minimum(c3),maximum(c3)))
  #display(p)

  c4 = reconstruction(bSF, b; lambd=0.01, SNRThresh = 3, frames=1:100, 
		      minFreq=80e3, nAverages=100, sparseTrafo="DCT", 
		      recChannels =1:2, iterations = 1, 
		      spectralLeakageCorrection=false, redFactor = redFactor,
		      useDFFoV = false)
  @test axisnames(c4) == names
  @test axisvalues(c4) == values
  p = imagesc(data(data(c4[1,:,:,1,1])), (minimum(c4),maximum(c4)))
  #display(p)

  c5 = reconstruction(bSF, b; lambd=0.01, SNRThresh = 3, frames=1:100,
		      minFreq=80e3, nAverages=100, sparseTrafo="DCT", 
		      recChannels =1:2, iterations = 1, 
		      spectralLeakageCorrection=true, redFactor = redFactor,
		      useDFFoV = true)
  @test axisnames(c5) == names
  @test axisvalues(c5) == valuesDF
  p = imagesc(data(data(c5[1,:,:,1,1])), (minimum(c5),maximum(c5)))
  #display(p)

  c6 = reconstruction(bSF, b; lambd=0.01, SNRThresh = 3, frames=1:100, 
		      minFreq=80e3, nAverages=100, sparseTrafo="DST", 
		      recChannels =1:2, iterations = 1, 		  
		      spectralLeakageCorrection=false, redFactor = redFactor, 
		      useDFFoV = false)
  @test axisnames(c6) == names
  @test axisvalues(c6) == values
  p = imagesc(data(data(c6[1,:,:,1,1])), (minimum(c6),maximum(c6)))
  #display(p)

  c7 = reconstruction(bSF, b; lambd=0.01, SNRThresh = 3, frames=1:100, 
		      minFreq=80e3, nAverages=100, sparseTrafo="DST", 
		      recChannels =1:2, iterations = 1, 
		      spectralLeakageCorrection=true, redFactor = redFactor,
		      useDFFoV = true)
  @test axisnames(c7) == names
  @test axisvalues(c7) == valuesDF
  p = imagesc(data(data(c7[1,:,:,1,1])), (minimum(c7),maximum(c7)))
  #display(p)
end

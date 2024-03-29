using MPIReco

@testset "sparse single- and multi-channel in-memory reconstruction" begin
  bSF = MPIFile(joinpath(datadir, "calibrations", "7.mdf"))
  b = MPIFile(joinpath(datadir, "measurements","20211226_204612_Dice", "1.mdf"))
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

  plan = getPlan("Single")
  setAll!(plan, :SNRThresh, 2)
  setAll!(plan, :frames, 1:100)
  setAll!(plan, :numAverages, 100)
  setAll!(plan, :minFreq, 80e3),
  setAll!(plan, :recChannels, 1:2)
  setAll!(plan, :iterations, 3)
  setAll!(plan, :spectralLeakageCorrection, false)
  setAll!(plan, :sf, bSF)
  setAll!(plan, :reg, [L2Regularization(0.1f0)])
  setAll!(plan, :solver, Kaczmarz)
  setAll!(plan, :gridding, SystemMatrixGriddingParameter(;gridsize=calibSize(bSF), fov = calibFov(bSF)))
  c1 = reconstruct(build(plan), b)
  @test axisnames(c1) == names
  @test axisvalues(c1) == values
  exportImage(joinpath(imgdir, "Sparse1.png"), Array(c1[1,:,:,1,1]))
  @test compareImg("Sparse1.png")

  plan = getPlan("Sparse")
  setAll!(plan, :SNRThresh, 2)
  setAll!(plan, :frames, 1:100)
  setAll!(plan, :minFreq, 80e3)
  setAll!(plan, :recChannels, 1:2)
  setAll!(plan, :iterations, 3)
  setAll!(plan, :spectralLeakageCorrection, false)
  setAll!(plan, :sf, bSF)
  setAll!(plan, :reg, [L2Regularization(0.1f0)])
  setAll!(plan, :numAverages, 100)
  setAll!(plan, :solver, Kaczmarz)
  #setAll!(plan, :tfCorrection, false)
  setAll!(plan, :redFactor, redFactor)
  
  setAll!(plan, :sparseTrafo, "FFT")
  setAll!(plan, :useDFFoV, false)
  c2 = reconstruct(build(plan), b)
  @test axisnames(c2) == names
  @test axisvalues(c2) == values
  exportImage(joinpath(imgdir, "Sparse2.png"), Array(c2[1,:,:,1,1]))
  @test compareImg("Sparse2.png")

  setAll!(plan, :sparseTrafo, "FFT")
  setAll!(plan, :useDFFoV, true)
  c3 = reconstruct(build(plan), b)
  @test axisnames(c3) == names
  @test axisvalues(c3) == valuesDF
  exportImage(joinpath(imgdir, "Sparse3.png"), Array(c3[1,:,:,1,1]))
  @test compareImg("Sparse3.png")

  setAll!(plan, :SNRThresh, 3)
  setAll!(plan, :iterations, 1)
  setAll!(plan, :reg, [L2Regularization(0.01f0)])
  setAll!(plan, :sparseTrafo, "DCT-IV")
  setAll!(plan, :useDFFoV, false)
  c4 = reconstruct(build(plan), b)
  @test axisnames(c4) == names
  @test axisvalues(c4) == values
  exportImage(joinpath(imgdir, "Sparse4.png"), Array(c4[1,:,:,1,1]))
  @test compareImg("Sparse4.png")

  setAll!(plan, :spectralLeakageCorrection, true)
  setAll!(plan, :sparseTrafo, "DCT-IV")
  setAll!(plan, :useDFFoV, true)
  c5 = reconstruct(build(plan), b)
  setAll!(plan, :spectralLeakageCorrection, false)
  @test axisnames(c5) == names
  @test axisvalues(c5) == valuesDF
  exportImage(joinpath(imgdir, "Sparse5.png"), Array(c5[1,:,:,1,1]))
  @test compareImg("Sparse5.png")

  setAll!(plan, :sparseTrafo, "DST")
  setAll!(plan, :useDFFoV, false)
  c6 = reconstruct(build(plan), b)
  @test axisnames(c6) == names
  @test axisvalues(c6) == values
  exportImage(joinpath(imgdir, "Sparse6.png"), Array(c6[1,:,:,1,1]))
  @test compareImg("Sparse6.png")

  setAll!(plan, :spectralLeakageCorrection, true)
  setAll!(plan, :sparseTrafo, "DST")
  setAll!(plan, :useDFFoV, true)
  c7 = reconstruct(build(plan), b)
  @test axisnames(c7) == names
  @test axisvalues(c7) == valuesDF
  exportImage(joinpath(imgdir, "Sparse7.png"), Array(c7[1,:,:,1,1]))
  @test compareImg("Sparse7.png")
end

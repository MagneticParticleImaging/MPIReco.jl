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

  params = Dict{Symbol, Any}()
  params[:SNRThresh] = 2
  params[:frames] = 1:100
  params[:numAverages] = 100
  params[:minFreq] = 80e3
  params[:recChannels] = 1:2
  params[:iterations] = 3
  params[:spectralLeakageCorrection] = false
  params[:sf] = bSF
  params[:reg] = [L2Regularization(0.1f0)]
  params[:solver] = Kaczmarz
  params[:gridding] = SystemMatrixGriddingParameter(;gridsize=calibSize(bSF), fov = calibFov(bSF))
  c1 = reconstruct("SinglePatch", b; params...)
  @test axisnames(c1) == names
  @test axisvalues(c1) == values
  exportImage(joinpath(imgdir, "Sparse1.png"), Array(c1[1,:,:,1,1]))
  @test compareImg("Sparse1.png")

  params = Dict{Symbol, Any}()
  params[:SNRThresh] = 2
  params[:frames] = 1:100
  params[:minFreq] = 80e3
  params[:recChannels] = 1:2
  params[:iterations] = 3
  params[:spectralLeakageCorrection] = false
  params[:sf] = bSF
  params[:reg] = [L2Regularization(0.1f0)]
  params[:numAverages] = 100
  params[:solver] = Kaczmarz
  #params[ :tfCorrection, false
  params[:redFactor] = redFactor

  c2 = reconstruct("SinglePatchSparse", b; params..., sparseTrafo = "FFT", useDFFoV = false)
  @test axisnames(c2) == names
  @test axisvalues(c2) == values
  exportImage(joinpath(imgdir, "Sparse2.png"), Array(c2[1,:,:,1,1]))
  @test compareImg("Sparse2.png")

  c3 = reconstruct("SinglePatchSparse", b; params..., sparseTrafo = "FFT", useDFFoV = true)
  @test axisnames(c3) == names
  @test axisvalues(c3) == valuesDF
  exportImage(joinpath(imgdir, "Sparse3.png"), Array(c3[1,:,:,1,1]))
  @test compareImg("Sparse3.png")

  params[:SNRThresh] = 3
  params[:iterations] = 1
  params[:reg] = [L2Regularization(0.01f0)]
  c4 = reconstruct("SinglePatchSparse", b; params..., sparseTrafo = "DCT-IV", useDFFoV = false)
  @test axisnames(c4) == names
  @test axisvalues(c4) == values
  exportImage(joinpath(imgdir, "Sparse4.png"), Array(c4[1,:,:,1,1]))
  @test compareImg("Sparse4.png")

  c5 = reconstruct("SinglePatchSparse", b; params..., sparseTrafo = "DCT-IV", useDFFoV = true, spectralLeakageCorrection = true)
  @test axisnames(c5) == names
  @test axisvalues(c5) == valuesDF
  exportImage(joinpath(imgdir, "Sparse5.png"), Array(c5[1,:,:,1,1]))
  @test compareImg("Sparse5.png")

  c6 = reconstruct("SinglePatchSparse", b; params..., sparseTrafo = "DST", useDFFoV = false)
  @test axisnames(c6) == names
  @test axisvalues(c6) == values
  exportImage(joinpath(imgdir, "Sparse6.png"), Array(c6[1,:,:,1,1]))
  @test compareImg("Sparse6.png")

  c7 = reconstruct("SinglePatchSparse", b; params..., sparseTrafo = "DST", useDFFoV = true, spectralLeakageCorrection = true)
  @test axisnames(c7) == names
  @test axisvalues(c7) == valuesDF
  exportImage(joinpath(imgdir, "Sparse7.png"), Array(c7[1,:,:,1,1]))
  @test compareImg("Sparse7.png")

  @testset "Multi Contrast" begin
    params = Dict{Symbol, Any}()
    params[:SNRThresh] = 2
    params[:frames] = 1:100
    params[:minFreq] = 80e3
    params[:recChannels] = 1:2
    params[:iterations] = 3
    params[:spectralLeakageCorrection] = false
    params[:sf] = MultiContrastFile([bSF, bSF])
    params[:reg] = [L2Regularization(0.1f0)]
    params[:numAverages] = 100
    params[:solver] = Kaczmarz
    params[:redFactor] = redFactor

    c = reconstruct("SinglePatchSparse", b; params..., sparseTrafo = "FFT", useDFFoV = false)
    @test size(c, 1) == length(params[:sf])
    @test size(c2)[2:end] == size(c)[2:end]

    c = reconstruct("SinglePatchSparse", b; params..., sparseTrafo = "FFT", useDFFoV = true)
    @test size(c, 1) == length(params[:sf])
    @test size(c3)[2:end] == size(c)[2:end]

    params[:SNRThresh] = 3
    params[:iterations] = 1
    params[:reg] = [L2Regularization(0.01f0)]
    c = reconstruct("SinglePatchSparse", b; params..., sparseTrafo = "DCT-IV", useDFFoV = false)
    @test size(c, 1) == length(params[:sf])
    @test size(c4)[2:end] == size(c)[2:end]

    c = reconstruct("SinglePatchSparse", b; params..., sparseTrafo = "DCT-IV", useDFFoV = true, spectralLeakageCorrection = true)
    @test size(c, 1) == length(params[:sf])
    @test size(c5)[2:end] == size(c)[2:end]

    c = reconstruct("SinglePatchSparse", b; params..., sparseTrafo = "DST", useDFFoV = false)
    @test size(c, 1) == length(params[:sf])
    @test size(c6)[2:end] == size(c)[2:end]

    c = reconstruct("SinglePatchSparse", b; params..., sparseTrafo = "DST", useDFFoV = true, spectralLeakageCorrection = true)
    @test size(c, 1) == length(params[:sf])
    @test size(c7)[2:end] == size(c)[2:end]

    @testset "Constructors" begin
      params = Dict{Symbol, Any}()
      params[:SNRThresh] = 2
      params[:frames] = 1:100
      params[:minFreq] = 80e3
      params[:recChannels] = 1:2
      params[:iterations] = 3
      params[:spectralLeakageCorrection] = false
      params[:sf] = MultiContrastFile([bSF, bSF])
      params[:reg] = [L2Regularization(0.1f0)]
      params[:numAverages] = 100
      params[:solver] = Kaczmarz
      params[:sparseTrafo] = "FFT"

      # Fitting scalar and vector values are accepted:
      plan = MPIRecoPlan("SinglePatchSparse"; params..., redFactor = 0.1)
      algo = build(plan)
      @test length(algo.S.parent.nzval)/prod(size(algo.S)) < 0.1

      plan = MPIRecoPlan("SinglePatchSparse"; params..., redFactor = [0.1, 0.1])
      algo = build(plan)
      @test length(algo.S.parent.nzval)/prod(size(algo.S)) < 0.1

      plan = MPIRecoPlan("SinglePatchSparse"; params..., redFactor = [0.1, 0.1, 0.1])
      @test_throws ArgumentError build(plan)

      plan = MPIRecoPlan("SinglePatchSparse"; params..., redFactor = [0.1, 0.05])
      algo = build(plan)
      @test length(algo.S.parent.nzval)/prod(size(algo.S)) < 0.08

      plan = MPIRecoPlan("SinglePatchSparse"; params..., redFactor = [0.05, 0.1])
      algo = build(plan)
      @test length(algo.S.parent.nzval)/prod(size(algo.S)) < 0.08

      # Same for threshold
      plan = MPIRecoPlan("SinglePatchSparse"; params..., thresh = 0.1)
      @test build(plan) isa SinglePatchReconstructionAlgorithm

      plan = MPIRecoPlan("SinglePatchSparse"; params..., thresh = [0.1, 0.1])
      @test build(plan) isa SinglePatchReconstructionAlgorithm

      plan = MPIRecoPlan("SinglePatchSparse"; params..., thresh = [0.1, 0.1, 0.1])
      @test_throws ArgumentError build(plan)

      plan = MPIRecoPlan("SinglePatchSparse"; params..., thresh = [0.1, 0.05])
      @test build(plan) isa SinglePatchReconstructionAlgorithm

      plan = MPIRecoPlan("SinglePatchSparse"; params..., thresh = [0.05, 0.1])
      @test build(plan) isa SinglePatchReconstructionAlgorithm
    end
  end
end

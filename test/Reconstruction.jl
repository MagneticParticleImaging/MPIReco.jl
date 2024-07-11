using MPIReco

@testset "single- and multi-channel in-memory reconstruction" begin
  bSF = MPIFile(joinpath(datadir, "calibrations", "12.mdf"))
  b = MPIFile(joinpath(datadir, "measurements", "20211226_203916_MultiPatch", "1.mdf"))
  names = (:color, :x, :y, :z, :time)
  values = (1:1,
	    -27.375u"mm":1.25u"mm":11.375u"mm",
	    -11.375u"mm":1.25u"mm":27.375u"mm",
	    0.0u"mm":1.0u"mm":0.0u"mm",
	    0.0u"ms":0.6528u"ms":0.0u"ms")

  # standard reconstruction
  params = Dict{Symbol, Any}()
  params[:SNRThresh] = 5
  params[:frames] = 1:1
  params[:minFreq] = 80e3
  params[:recChannels] = 1:2
  params[:iterations] = 1
  params[:spectralLeakageCorrection] = true
  params[:sf] = bSF
  params[:reg] = [L2Regularization(0.0f0)]
  params[:solver] = Kaczmarz

  c1 = reconstruct("SinglePatch", b; params...)
  @test axisnames(c1) == names
  @test axisvalues(c1) == values
  exportImage(joinpath(imgdir, "Reconstruction1.png"), Array(c1[1,:,:,1,1]))
  @test compareImg("Reconstruction1.png")

  # fused lasso
  #params[:iterations, 100)
  #params[:loadasreal, true)
  #params[:solver, FusedLasso)
  #params[:reg, [TVRegularization(0.01f0), L1Regularization(0.01f0)])
  #params[:normalizeReg, NoNormalization())
  #c2 = reconstruct(build(plan), b)
  #params[:iterations, 1)
  #params[:loadasreal, false)
  #params[:solver, Kaczmarz)
  #params[:reg, [L2Regularization(0.0f0)])
  #params[:normalizeReg, SystemMatrixBasedNormalization())
  #@test axisnames(c2) == names
  #@test axisvalues(c2) == values
  #exportImage(joinpath(imgdir, "Reconstruction2.png"), Array(c2[1,:,:,1,1]))
  #@test compareImg("Reconstruction2.png")

  # with interpolation
  params[:gridding] = SystemMatrixGriddingParameter(;gridsize=[100,100,1], fov = calibFov(bSF))
  c3 = reconstruct("SinglePatch", b; params...)
  @test axisnames(c3) == names
  @test axisvalues(c3) == (values[1], -27.8u"mm":0.4u"mm":11.8u"mm", -11.8u"mm":0.4u"mm":27.8u"mm", values[4:5]...)
  exportImage(joinpath(imgdir, "Reconstruction3.png"), Array(c3[1,:,:,1,1]))
  @test compareImg("Reconstruction3.png")

  # with fov adpation and center shift
  params[:gridding] = SystemMatrixGriddingParameter(;gridsize=[50,50,1], fov = [0.04,0.04,1], center = [0.0,-0.01,0])
  c4 = reconstruct("SinglePatch", b; params...)
  @test axisnames(c4) == names
  # TODO Tobi: does this make sense?
  @test axisvalues(c4) == (values[1], -27.6u"mm":0.8u"mm":11.6u"mm", -11.6u"mm":0.8u"mm":27.6u"mm", 499.5u"mm":1000.0u"mm":499.5u"mm", values[5])
  exportImage(joinpath(imgdir, "Reconstruction4.png"), Array(c4[1,:,:,1,1]))
  @test compareImg("Reconstruction4.png")

  params[:gridding] = SystemMatrixGriddingParameter(;gridsize=calibSize(bSF), fov = calibFov(bSF))

  # multi colored reconstruction (new interface)
  params[:sf] = MultiContrastFile([bSF,bSF])
  c5 = reconstruct("SinglePatch", b; params...)
  @test axisnames(c5) == names
  @test axisvalues(c5) == (1:2,values[2:end]...)
  exportImage(joinpath(imgdir, "Reconstruction5a.png"), Array(c5[1,:,:,1,1]))
  @test compareImg("Reconstruction5a.png")
  exportImage(joinpath(imgdir, "Reconstruction5b.png"), Array(c5[2,:,:,1,1]))
  @test compareImg("Reconstruction5b.png")
#=
  # dict based reco
  r = defaultRecoParams()
  r[:measPath] = filepath(b)
  r[:SFPath] = filepath(bSF)
  r[:frames] =  1
  r[:minFreq] = 80e3
  r[:SNRThresh] = 4
  r[:lambd] = 0.001
  r[:iterations] = 1
  c6 = reconstruction(r)
  @test axisnames(c6) == names
  @test axisvalues(c6) == values
  exportImage(joinpath(imgdir, "Reconstruction6.png"), Array(c6[1,:,:,1,1]))
  @test compareImg("Reconstruction6.png")
  =#

  params[:sf] = bSF
  params[:reg] = [L2Regularization(0.1f0)]
  params[:gridding] = SystemMatrixGriddingParameter(;gridsize=calibSize(bSF), fov = calibFov(bSF))
  # channel weighting
  params[:weightingParams] = ChannelWeightingParameters(channelWeights = [1.0, 1.0, 1.0])
  c7a = reconstruct("SinglePatch", b; params...)
  exportImage(joinpath(imgdir, "Reconstruction7a.png"), Array(c7a[1,:,:,1,1]))
  @test compareImg("Reconstruction7a.png")

  params[:weightingParams] = ChannelWeightingParameters(channelWeights = [1.0,0.001,1.0])
  c7b = reconstruct("SinglePatch", b; params...)
  exportImage(joinpath(imgdir, "Reconstruction7b.png"), Array(c7b[1,:,:,1,1]))
  @test compareImg("Reconstruction7b.png")

  params[:weightingParams] = ChannelWeightingParameters(channelWeights = [0.001,1.0,1.0])
  c7c = reconstruct("SinglePatch", b; params...)
  exportImage(joinpath(imgdir, "Reconstruction7c.png"), Array(c7c[1,:,:,1,1]))
  @test compareImg("Reconstruction7c.png")
  

  params[:weightingParams] = ChannelWeightingParameters(channelWeights = [1.0, 0.0, 1.0])
  c7d = reconstruct("SinglePatch", b; params...)
  params[:weightingParams] = NoWeightingParameters()
  params[:recChannel] = 1:1
  c7e = reconstruct("SinglePatch", b; params...)
  params[:recChannel] = 1:2
  @test isapprox(arraydata(c7d), arraydata(c7e))

  params[:weightingParams] = ChannelWeightingParameters(channelWeights = [0.0, 1.0, 1.0])
  c7f = reconstruct("SinglePatch", b; params...)
  params[:weightingParams] = NoWeightingParameters()
  params[:recChannel] = 2:2
  c7g = reconstruct("SinglePatch", b; params...)
  params[:recChannel] = 1:2
  @test isapprox(arraydata(c7f), arraydata(c7g))

  params[:weightingParams] = WhiteningWeightingParameters(whiteningMeas = bSF)
  c8 = reconstruct("SinglePatch", b; params...)
  exportImage(joinpath(imgdir, "Reconstruction8.png"), Array(c8[1,:,:,1,1]))
  @test compareImg("Reconstruction8.png")
end

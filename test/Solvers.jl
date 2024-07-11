using MPIReco

@testset "Different Solver Reconstructions" begin
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
  params[:spectralLeakageCorrection] = true
  params[:sf] = bSF
  params[:gridding] = SystemMatrixGriddingParameter(;gridsize=calibSize(bSF), fov = calibFov(bSF))
  params[:weightingParams] = WhiteningWeightingParameters(whiteningMeas = bSF)


  params[:reg] = [L2Regularization(0.1f0)]
  params[:iterations] = 100

  for solver in [Kaczmarz, CGNR, FISTA, OptISTA, POGM] # , ADMM, SplitBregman]
    params[:solver] = solver
    params[:rho] = 0.3f0
    c = reconstruct("SinglePatch", b; params...)
    @test axisnames(c) == names
    @test axisvalues(c) == values
    exportImage(joinpath(imgdir, "Solver_$solver.png"), Array(c[1,:,:,1,1]))
    @test compareImg("Solver_$solver.png")
  end

end

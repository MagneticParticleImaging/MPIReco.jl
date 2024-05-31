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
  plan = getPlan("Single")
  plan.parameter.reco.solverParams = RecoPlan(ElaborateSolverParameters)
  setAll!(plan, :SNRThresh, 5)
  setAll!(plan, :frames, 1:1)
  setAll!(plan, :minFreq, 80e3),
  setAll!(plan, :recChannels, 1:2)
  setAll!(plan, :spectralLeakageCorrection, true)
  setAll!(plan, :sf, bSF)
  setAll!(plan, :gridding, SystemMatrixGriddingParameter(;gridsize=calibSize(bSF), fov = calibFov(bSF)))
  setAll!(plan, :weightingParams, WhiteningWeightingParameters(whiteningMeas = bSF))


  setAll!(plan, :reg, [L2Regularization(0.1f0)])
  setAll!(plan, :iterations, 100)

  for solver in [Kaczmarz, CGNR, FISTA, OptISTA, POGM] # , ADMM, SplitBregman]
    setAll!(plan, :solver, solver)
    setAll!(plan, :rho, 0.3f0)
    c = reconstruct(build(plan), b)
    @test axisnames(c) == names
    @test axisvalues(c) == values
    exportImage(joinpath(imgdir, "Solver_$solver.png"), Array(c[1,:,:,1,1]))
    @test compareImg("Solver_$solver.png")
  end

end

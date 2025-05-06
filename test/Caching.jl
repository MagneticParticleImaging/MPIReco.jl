@testset "Plan Caching" begin
  emptyRecoCache!()
  @test isempty(MPIReco.recoPlans)

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

  # Default is no caching
  reconstruct("SinglePatch", b; params...)
  @test isempty(MPIReco.recoPlans)

  # Enable caching
  reconstruct("SinglePatch", b, true; params...)
  @test !isempty(MPIReco.recoPlans)

  # Cache clearing (Kinda hacky because the testset depends on this already)
  emptyRecoCache!()
  @test isempty(MPIReco.recoPlans)

  # Distributed
  db = DMPIFile(joinpath(datadir, "measurements", "20211226_203916_MultiPatch", "1.mdf"), worker = 1)
  reconstruct("SinglePatch", db; params...)
  @test isempty(MPIReco.recoPlans)

  reconstruct("SinglePatch", db, true; params...)
  @test !isempty(MPIReco.recoPlans)

end
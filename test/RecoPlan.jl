@testset "RecoPlan" begin
  @testset "Construction" begin
    bSF = MPIFile(joinpath(datadir, "calibrations", "12.mdf"))
    b = MPIFile(joinpath(datadir, "measurements", "20211226_203916_MultiPatch", "1.mdf"))
    
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

    plan = MPIRecoPlan("SinglePatch"; params...)
    algo = build(plan)
    @test algo isa SinglePatchReconstructionAlgorithm
    @test isapprox(arraydata(reconstruct(algo, b)), arraydata(reconstruct("SinglePatch", b; params...)))

    setAll!(plan, :weightingParams, WhiteningWeightingParameters(bSF, false))
    algo = build(plan)
    reconstruct(algo, b)
    @test !isnothing(algo.weights)
  end
end
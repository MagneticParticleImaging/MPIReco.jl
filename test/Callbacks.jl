@testset "Callbacks" begin

  bSF = MPIFile(joinpath(datadir, "calibrations", "12.mdf"))
  b = MPIFile(joinpath(datadir, "measurements", "20211226_203916_MultiPatch", "1.mdf"))
  params = Dict{Symbol, Any}()
  params[:SNRThresh] = 5
  params[:frames] = 1:1
  params[:minFreq] = 80e3
  params[:recChannels] = 1:2
  params[:iterations] = 2
  params[:spectralLeakageCorrection] = true
  params[:sf] = bSF
  params[:reg] = [L2Regularization(0.0f0)]
  params[:solver] = Kaczmarz
  
    
  reconstruct("SinglePatch", b; params...) do solver, frame, iteration
    @test solver isa params[:solver]
    @test frame == 1
  end

  params[:frames] = [1, 1, 1]
  cb = StoreSolutionPerFrameCallback()
  c = reconstruct(cb, "SinglePatch", b; params...)
  @test length(cb.solutions) == length(params[:frames])
  @test map(length, cb.solutions) == fill(params[:iterations], 3) .+ 1
  @test all(map(solutions -> iszero(solutions[1]), cb.solutions))
  @test all(map(solutions -> isapprox(solutions[end], vec(c[1, :, :, :, 1])), cb.solutions))

  ref = reshape(c[1, :, :, :, :].data.data, :, 3)
  cb = CompareSolutionPerFrameCallback(ref)
  c2 = reconstruct(cb, "SinglePatch", b; params...)
  @test length(cb.results) == length(params[:frames])
end
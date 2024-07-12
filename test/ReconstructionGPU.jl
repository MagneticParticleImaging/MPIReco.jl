for arrayType in arrayTypes
  @testset "Reconstruction $arrayType" begin
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
    params[:solver] = CGNR
  
    c1 = reconstruct("SinglePatch", b; params...)
    @test axisnames(c1) == names
    @test axisvalues(c1) == values
  
    c2 = reconstruct("SinglePatch", b; params..., arrayType = arrayType)
    @test axisnames(c2) == names
    @test axisvalues(c2) == values

    @test isapprox(arraydata(c1), arraydata(c2))
  end
end
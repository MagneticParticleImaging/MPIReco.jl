using MPIReco

@testset "SystemMatrixExtrapolation" begin
    bSF = MPIFile(joinpath(datadir, "calibrations", "12.mdf"))
    freq = sort(filterFrequencies(bSF,SNRThresh=5,minFreq=60e3,recChannels=1:2))
    SMextr1 = extrapolateSM(bSF, freq, (5,5,0); method=1)
    
    SM, grid = getSF(bSF, freq, returnasmatrix=true)
    SMextr2 = extrapolateSM(SM, grid, (5,5,0); method=1)
    SMextr3 = extrapolateSM(SM, grid, 5)
    
    SMsizeExtr = shape(grid) .+ [10;10;0]
    @test SMextr2 == transpose(SMextr1)
    @test SMextr2 == SMextr3
    @test prod(SMsizeExtr) == size(SMextr2,1)
end
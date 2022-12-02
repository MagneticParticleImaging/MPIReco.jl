using MPIReco

@testset "SystemMatrixExtrapolation" begin
    bSF = MPIFile(joinpath(datadir, "calibrations", "12.mdf"))
    freq = sort(filterFrequencies(bSF,SNRThresh=5,minFreq=60e3,recChannels=1:2))
    SMextr1 = extrapolateSM(bSF, freq, (5,5,0); method=1)[1]
    
    SM, grid = getSF(bSF, freq, returnasmatrix=true)
    SMextr2, extrgrid = extrapolateSM(SM, grid, (5,5,0); method=1)
    SMextr3 = extrapolateSM(SM, grid, 5)[1]
    
    @test SMextr2 == transpose(SMextr1)
    @test SMextr2 == SMextr3

    SMsizeExtr = shape(grid) .+ [10;10;0]
    @test shape(extrgrid) == SMsizeExtr
    @test prod(SMsizeExtr) == size(SMextr2,1)
    
    b = MPIFile(joinpath(datadir, "measurements", "20211226_203916_MultiPatch", "1.mdf"))
    c1 = reconstruction(bSF, b; SNRThresh=5, frames=1, minFreq=80e3,
    recChannels=1:2, iterations=1, spectralLeakageCorrection=true, SMextrapolation=nothing)    
    c_extr = reconstruction(bSF, b; SNRThresh=5, frames=1, minFreq=80e3,
    recChannels=1:2, iterations=1, spectralLeakageCorrection=true, SMextrapolation=(3,3))
    @test size(c1[1,:,:,:,1]) .+ (6,6,0) == size(c_extr[1,:,:,:,1]) 
    exportImage(joinpath(imgdir, "Extrapolated1.png"), arraydata(c_extr[1,:,:,1,1]))

    SFdirs = ["8.mdf", "9.mdf", "10.mdf", "11.mdf"]
    bSFs = MultiMPIFile(joinpath.(datadir, "calibrations", SFdirs))
    bdirs = ["1.mdf", "2.mdf", "3.mdf", "4.mdf"]
    b = MultiMPIFile(joinpath.(datadir, "measurements", "20211226_203916_MultiPatch", bdirs))
    c2 = reconstruction(bSFs, b;
                  SNRThresh=5, frames=1, minFreq=80e3,
                  recChannels=1:2,iterations=1,
                  spectralLeakageCorrection=false
                  )
    c_extr2 = reconstruction(bSFs, b;
                  SNRThresh=5, frames=1, minFreq=80e3,
                  recChannels=1:2,iterations=1,
                  spectralLeakageCorrection=false,
                  SMextrapolation=3)                  
    @test size(c2[1,:,:,:,1]) .+ (6,6,0) == size(c_extr2[1,:,:,:,1])               
    exportImage(joinpath(imgdir, "ExtrapolatedMultiPatch1.png"), arraydata(c_extr2[1,:,:,1,1]))
end
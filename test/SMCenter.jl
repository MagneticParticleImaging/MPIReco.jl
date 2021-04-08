using MPIReco

@testset "system matrix center estimation" begin
  SFs = ["$datadir/SF_MP01", "$datadir/SF_MP02", "$datadir/SF_MP03", "$datadir/SF_MP04"]
  center = [[9.0,23.0,1.0],[9.0,10.0,1.0],[22.0,23.0,1.0],[22.0,10.0,1.0]]

  for (l,SF) in enumerate(SFs)
    bSF = MPIFile(SF)
    a = findCenterOfDfFov(bSF)
    @test a == center[l]
    S = getSF(bSF,2,0,0,2)
    exportImage("$imgdir/Center1.png", abs.(S)[:,:,1])
    #p = imagesc(abs.(S)[:,:,1])
    #add(p, Points([a[2]-1],[a[1]-1],color="yellow",lw=5))
    #savefig(p, "$imgdir/Center1.png")
  end

  for (l,SF) in enumerate(SFs)
    bSF = MPIFile(SF)
    a = findCenterOfDfFov(bSF)
    @test a == center[l]
    S = getSF(bSF,5,6,0,1)
    exportImage("$imgdir/Center2.png", abs.(S)[:,:,1])
    #p = imagesc(abs.(S)[:,:,1])
    #add(p, Points([a[2]-1],[a[1]-1],color="yellow",lw=5))
    #savefig(p, "$imgdir/Center2.png")
  end
end

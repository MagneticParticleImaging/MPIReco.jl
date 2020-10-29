using MPIReco
using Test

@testset "test saving and loading of MPI images to and from MDF's" begin
  bSF = MPIFile("$datadir/SF_MP")
  b = MPIFile("$datadir/MP01")
  r = defaultRecoParams()
  r[:measPath] = filepath(b)
  r[:SFPath] = filepath(bSF)
  r[:frames] =  1
  r[:minFreq] = 80e3
  r[:SNRThresh] = 4
  r[:lambd] = 0.001
  r[:iterations] = 1

  # standard reconstruction
  c = reconstruction(r)
  names =  axisnames(c)
  values =  axisvalues(c)

  # save c as MDF
  saveRecoData("$imgdir/reco.mdf",c)

  # load MDF
  cmdf = loadRecoData("$imgdir/reco.mdf")
  @test axisnames(cmdf) == names
  @test axisvalues(cmdf) == values
end

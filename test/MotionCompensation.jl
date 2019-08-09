using MPIReco

@testset "motion compensation reconstruction" begin

  bSF = MultiMPIFile(["motionComp/SF1Small.mdf", "motionComp/SF2Small.mdf",
                      "motionComp/SF3Small.mdf", "motionComp/SF4Small.mdf"])
  bBG = MPIFile("motionComp/measBG.mdf")
  b = MPIFile("motionComp/measFast.mdf") # "measSlow.mdf"

end

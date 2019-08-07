using MPIReco

@testset "motion compensation reconstruction" begin

  bSF = MultiMPIFile(["SF1Small.mdf", "SF2Small.mdf", "SF3Small.mdf", "SF4Small.mdf"])
  bBG = MPIFile("measBG.mdf")
  b = MPIFile("measFast.mdf") # "measSlow.mdf"

end

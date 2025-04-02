using MPIReco

@testset "Cartesian reconstruction" begin
  fnSM = joinpath(datadir, "calibrations", "1.mdf")
  fnSMProc1 = joinpath(datadir, "calibrations", "1p1.mdf")
  bSF = MPIFile(fnSM)
  b = MPIFile(joinpath(datadir, "measurements", "20211226_203850_HeadScanner", "1.mdf"))

  numPeriodAverages = 65
  numPatches = div(acqNumPeriodsPerFrame(bSF), numPeriodAverages)
  bgCorrection = false

  # MP Reco # TODO this does not seem to branch into mp reco atm
  numPeriodGrouping = 1
  maxMixingOrder = -1

  plan = loadPlan(joinpath(@__DIR__(), "..", "config", "SinglePatch"), [MPIReco, RegularizedLeastSquares, MPIFiles, AbstractImageReconstruction])
  setAll!(plan, :sf, bSF)
  setAll!(plan, :frames, 1:10)
  setAll!(plan, :numAverages, 10)
  setAll!(plan, :numPeriodGrouping, numPeriodGrouping)
  setAll!(plan, :numPeriodAverages, numPeriodAverages)
  setAll!(plan, :minFreq, 30e3)
  setAll!(plan, :maxFreq, 600e3)
  #setAll!(plan, :maxMixingOrder, maxMixingOrder) # This was used in old tests, but didn't have any effect anymore, now it does
  setAll!(plan, :solver, Kaczmarz)
  setAll!(plan, :reg, [L2Regularization(0.01f0)])
  setAll!(plan, :iterations, 100)
  setAll!(plan, :bgParams, NoBackgroundCorrectionParameters())
  setAll!(plan, :spectralLeakageCorrection, false)
  setAll!(plan, :gridsize, calibSize(bSF))
  setAll!(plan, :fov, calibFov(bSF))  
  @time c1 = reconstruct(build(plan), b)

  exportImage(joinpath(imgdir, "Cartesian1.png"), Array(c1[1,:,:,1,1]))
  @test compareImg("Cartesian1.png")

  # SP Reco
  setAll!(plan, :numPeriodGrouping, numPatches)
  #setAll!(plan, :maxMixingOrder, 12) # see above
  @time c2 = reconstruct(build(plan), b)
  exportImage(joinpath(imgdir, "Cartesian2.png"), Array(c2[1,:,:,1,1]))
  @test compareImg("Cartesian2.png") skip = true

  # from Postprocessed
  saveasMDF(fnSMProc1, fnSM, numPeriodAverages=65, applyCalibPostprocessing=true, numPeriodGrouping=100)
  bSFProc1 = MPIFile(fnSMProc1)

  setAll!(plan, :sf, bSFProc1)
  setAll!(plan, :numPeriodGrouping, 1)
  setAll!(plan, :numPeriodAverages, 1)
  setAll!(plan.parameter.pre, :numPeriodGrouping, 100) # SM is preprocessed, only need to do it for meas data
  setAll!(plan.parameter.pre, :numPeriodAverages, 65)
  @time c3 = reconstruct(build(plan), b)

  exportImage(joinpath(imgdir, "Cartesian3.png"), Array(c3[1,:,:,1,1]))
  @test compareImg("Cartesian3.png") skip = true

  ####  Low Level ####

  N = calibSize(bSF)

  ## MP
  numPeriodGrouping = 1
  maxMixingOrder = -1

  freq = filterFrequencies(bSF, numPeriodGrouping=numPeriodGrouping, numPeriodAverages=numPeriodAverages,
             minFreq=30e3, maxFreq=600e3,  maxMixingOrder=maxMixingOrder)

  S, grid = getSF(bSF, freq, nothing, "Kaczmarz", numPeriodAverages=numPeriodAverages, bgCorrection=bgCorrection,
           numPeriodGrouping=numPeriodGrouping, spectralLeakageCorrection=false);

  u = getMeasurementsFD(b, frequencies=freq, bgCorrection=bgCorrection, frames=1:10, numAverages=10,
               numPeriodGrouping=numPeriodGrouping, numPeriodAverages=numPeriodAverages,
               spectralLeakageCorrection=false)

  @time c4 = reshape(reconstruction(S, u, reg = [L2Regularization(0.01f0), PositiveRegularization()], iterations=100), N[1], N[2])

  exportImage(joinpath(imgdir, "Cartesian4.png"), c4)
  @test compareImg("Cartesian4.png") skip = true

  ## SP
  numPeriodGrouping = numPatches
  maxMixingOrder = 12

  freq = filterFrequencies(bSF, numPeriodGrouping=numPeriodGrouping, numPeriodAverages=numPeriodAverages,
             minFreq=30e3, maxFreq=600e3, numSidebandFreqs = 12)

  S, grid = getSF(bSF, freq, nothing, "Kaczmarz", numPeriodAverages=numPeriodAverages, bgCorrection=bgCorrection,
           numPeriodGrouping=numPeriodGrouping, spectralLeakageCorrection=false);

  u = getMeasurementsFD(b, frequencies=freq, bgCorrection=bgCorrection, frames=1:10, numAverages=10,
               numPeriodGrouping=numPeriodGrouping, numPeriodAverages=numPeriodAverages,
               spectralLeakageCorrection=false)

  @time c5 = reshape(reconstruction(S, u, reg = [L2Regularization(0.01f0), PositiveRegularization()], iterations=100), N[1], N[2])

  exportImage(joinpath(imgdir, "Cartesian5.png"), c5)
  @test compareImg("Cartesian5.png") skip = true

  ####  Multi Color ####

  setAll!(plan, :sf, MultiContrastFile([bSF,bSF]))
  setAll!(plan, :numPeriodGrouping, numPeriodGrouping)
  setAll!(plan, :numPeriodAverages, numPeriodAverages)
  @time c6 = reconstruct(build(plan), b)

  exportImage(joinpath(imgdir, "Cartesian6.png"), Array(c6[1,:,:,1,1]))
  @test compareImg("Cartesian6.png") skip = true

end

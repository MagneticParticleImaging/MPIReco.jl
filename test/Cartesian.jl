using MPIReco

@testset "sparse single- and multi-channel in-memory reconstruction" begin
  fnSM = joinpath(datadir, "calibrations", "1.mdf")
  fnSMProc1 = joinpath(datadir, "calibrations", "1p1.mdf")
  bSF = MPIFile(fnSM)
  b = MPIFile(joinpath(datadir, "measurements", "20211226_203850_HeadScanner", "1.mdf"))

  numPeriodAverages = 65
  numPatches = div(acqNumPeriodsPerFrame(bSF), numPeriodAverages)
  bgCorrection = false

  # MP Reco
  numPeriodGrouping = 1
  maxMixingOrder = -1

  @time c1 = reconstruction(bSF, b, frames=1:10, numAverages=10,
           numPeriodGrouping=numPeriodGrouping, numPeriodAverages=numPeriodAverages,
           minFreq=30e3, maxFreq=600e3, maxMixingOrder=maxMixingOrder,
           lambd=0.01, iterations=100, bgCorrectionInternal=bgCorrection,
           spectralLeakageCorrection=false)

  exportImage(joinpath(imgdir, "Cartesian1.png"), arraydata(c1[1,:,:,1,1]))

  # SP Reco
  numPeriodGrouping = numPatches
  maxMixingOrder = 12

  @time c2 = reconstruction(bSF, b, frames=1:10, numAverages=10,
           numPeriodGrouping=numPeriodGrouping, numPeriodAverages=numPeriodAverages,
           minFreq=30e3, maxFreq=600e3, maxMixingOrder=maxMixingOrder,
           lambd=0.01, iterations=100, bgCorrectionInternal=bgCorrection,
           spectralLeakageCorrection=false)

  exportImage(joinpath(imgdir, "Cartesian2.png"), arraydata(c2[1,:,:,1,1]))

  # from Postprocessed
  saveasMDF(fnSMProc1, fnSM, numPeriodAverages=65, applyCalibPostprocessing=true, numPeriodGrouping=100)
  bSFProc1 = MPIFile(fnSMProc1)

  maxMixingOrder = 12

  @time c2 = reconstruction(bSFProc1, b, frames=1:10, numAverages=10,
           numPeriodGrouping=1, numPeriodAverages=1, SNRThresh=0,
           minFreq=30e3, maxFreq=600e3, #maxMixingOrder=maxMixingOrder,
           lambd=0.01, iterations=100, bgCorrectionInternal=bgCorrection,
           spectralLeakageCorrection=false)

  exportImage(joinpath(imgdir, "Cartesian3.png"), arraydata(c2[1,:,:,1,1]))

  ####  Low Level ####

  N = calibSize(bSF)

  ## MP
  numPeriodGrouping = 1
  maxMixingOrder = -1

  freq = filterFrequencies(bSF, numPeriodGrouping=numPeriodGrouping, numPeriodAverages=numPeriodAverages,
             minFreq=30e3, maxFreq=600e3,  maxMixingOrder=maxMixingOrder)

  S, grid = getSF(bSF, freq, numPeriodAverages=numPeriodAverages, bgCorrection=bgCorrection,
           numPeriodGrouping=numPeriodGrouping, spectralLeakageCorrection=false);

  u = getMeasurementsFD(b, frequencies=freq, bgCorrection=bgCorrection, frames=1:10, numAverages=10,
               numPeriodGrouping=numPeriodGrouping, numPeriodAverages=numPeriodAverages,
               spectralLeakageCorrection=false)

  c3 = reshape(reconstruction(transpose(S), u, lambd=0.01, iterations=100), N[1], N[2])

  exportImage(joinpath(imgdir, "Cartesian4.png"), c3)


  ## SP
  numPeriodGrouping = numPatches
  maxMixingOrder = 12

  freq = filterFrequencies(bSF, numPeriodGrouping=numPeriodGrouping, numPeriodAverages=numPeriodAverages,
             minFreq=30e3, maxFreq=600e3,  maxMixingOrder=maxMixingOrder)

  S, grid = getSF(bSF, freq, numPeriodAverages=numPeriodAverages, bgCorrection=bgCorrection,
           numPeriodGrouping=numPeriodGrouping, spectralLeakageCorrection=false);

  u = getMeasurementsFD(b, frequencies=freq, bgCorrection=bgCorrection, frames=1:10, numAverages=10,
               numPeriodGrouping=numPeriodGrouping, numPeriodAverages=numPeriodAverages,
               spectralLeakageCorrection=false)

  c4 = reshape(reconstruction(transpose(S), u, lambd=0.01, iterations=100), N[1], N[2])

  exportImage(joinpath(imgdir, "Cartesian5.png"), c4)

  ####  Multi Color ####

  @time c5 = reconstruction(MultiContrastFile([bSF,bSF]), b, frames=1:10, numAverages=10,
           numPeriodGrouping=numPeriodGrouping, numPeriodAverages=numPeriodAverages,
           minFreq=30e3, maxFreq=600e3, maxMixingOrder=maxMixingOrder,
           lambd=0.01, iterations=100, bgCorrectionInternal=bgCorrection,
           spectralLeakageCorrection=false)

  exportImage(joinpath(imgdir, "Cartesian6.png"), arraydata(c5[1,:,:,1,1]))
  
end
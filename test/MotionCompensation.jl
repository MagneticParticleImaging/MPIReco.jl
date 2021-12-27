using MPIReco

@testset "motion compensation reconstruction" begin

  ################
  ## Parameters ##
  ################
  alpha = 1.2

  # Choose the peak number according to the expected approximate value. Be careful, if there are peaks not corresponding to object motion, division by choosePeak can be false.
  # For experiments with different periodic motions it is used to pick one motion frequency (then don't divide by choosePeak)

  choosePeak = 4

  windowType = 1  # 1: Hann, 2: FT1A05, 3: Rectangle

  # Measurement data
  datadirMeas = joinpath(datadir, "measurements", "20211226_204336_Rotation")
  bMeas = MPIFile(joinpath(datadirMeas, "1.mdf")) # high frequency
  bBG = MPIFile(joinpath(datadirMeas, "2.mdf")) # background measurement

  # System matrices
  datadirSF = joinpath(datadir, "calibrations")
  SFall = ["3.mdf","4.mdf","5.mdf","6.mdf"]
  bSF = MultiMPIFile(joinpath.(datadirSF, SFall))

  # Background frames
  bgFrames = [1:200, 201:400, 401:600 ,601:800]

  # Selected frames
  recoFrame = 1

  ####################
  ## Reconstruction ##
  ####################

  # Reconstruction parameters
  lambda = 0.01
  iterations = 1

  c = reconstruction(bSF, bMeas, periodicMotionCorrection = true,
                                bEmpty = bBG, bgFrames = bgFrames, #background measurement
                                alpha = alpha, choosePeak = choosePeak, recoFrame = recoFrame,
                                windowType=windowType, higherHarmonic=choosePeak,
                                lambda=lambda, iterations=iterations, # reconstruction parameter
                                SNRThresh=10, minFreq=80e3, recChannels=[1,2,3]) #frequency selection

  exportImage(joinpath(imgdir, "MotionComp.png"), maximum(arraydata(c[1,:,:,:,1]),dims=3)[:,:,1])

end

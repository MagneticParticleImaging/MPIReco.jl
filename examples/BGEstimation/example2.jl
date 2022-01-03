using MPIReco, ImageView #, PyPlot

# Download data
# include("downloadData.jl")

# Measurement data
datadir = "/Users/knopp/.julia/artifacts/dc0fe8f812428df92eae484f90a9d2d2c91261fe/"
datadirMeas = joinpath(datadir, "./measurements/20211226_204133_Background/")
f = MPIFile(datadirMeas*"1.mdf")
fEmpty = MPIFile(datadirMeas*"2.mdf")

#fgFrames = vcat([1], collect(17:35), collect(49:69), collect(83:103), 
#                collect(117:137), collect(151:171) ,[200])
            
fgFrames = vcat([1,17,18,19,33,34,35,49,50,51,67,68,69,83,84,85,101,102,103,
                 115,116,117, 151,152,153, 169,170,171, 200])

# System matrices
datadirSF = joinpath(datadir, "./calibrations/")
fSF = MPIFile(datadirSF*"2.mdf")

# Put all parameters into a dictionary
params = Dict{Symbol,Any}(
  :frames => 1:200,
  :emptyMeas => fEmpty,
  #:bgFrames => 1:20,
  #:emptyMeas => f,
  :bgFrames => 1:19,
  :numAverages => 1,
  :bgCorrectionInternal => false,
  :λ => 1e-2 / 1, # 1e-1
  :recChannels => [1,2],
  :iterations => 20,
  :minFreq => 80e3,
  :SNRThresh => 3.0,
  :enforceReal => true,
  :enforcePositive => true,
  :spectralLeakageCorrection =>false
 )

############################
## Perform reconstruction ##
############################

@info "Reconstruction with static BG correction"
paramsStatic = deepcopy(params)
#cStaticCorr = reconstruction(fSF,f; paramsStatic...).data.data


@info "Reconstruction with joint estimation of BG and particle distribution"
paramsJoint = deepcopy(params)
paramsJoint[:bgDictSize] = 10
paramsJoint[:β] = 2.56e-7
paramsJoint[:λ] = 1e-2 
paramsJoint[:reco] = :tempReg
paramsJoint[:numBGAverages] = 1
cBGEst = reconstruction(fSF,f; paramsJoint...).data.data

@info "Reconstruction with joint estimation of BG and particle distribution"
paramsTRTikhonov = deepcopy(params)
paramsTRTikhonov[:bgDictSize] = 10
paramsTRTikhonov[:β] = [2.56e-7, 2.56e-7] 
paramsTRTikhonov[:λ] = [1e-2, 1e-2]
paramsTRTikhonov[:reco] = :tempReg
paramsTRTikhonov[:numBGAverages] = 1
cTRTikhonov = reconstruction(fSF,f; paramsTRTikhonov...).data.data

#ImageView.imshow(abs.(cat(cBGEst, cTRTikhonov, cTRTikhonov.- cBGEst,  dims=3)[1,:,:,1,:]),CLim(0,maximum(cBGEst)), name="TempReg")


L = 200

@info "Temp Regu Reco"
paramsTempReg = deepcopy(params)
paramsTempReg[:reco] = :tempReg
paramsTempReg[:bgDictSize] = 10
paramsTempReg[:β] = 2.56e-7 #/ 1000
paramsTempReg[:λ] = 1e-2  #* 100 # / 1000
paramsTempReg[:numBGAverages] = 1
#paramsTempReg[:idxFG] = fgFrames
paramsTempReg[:idxFG] = push!(collect(1:2:(L-1)),L) # [1, 2, 5, 9,  10] #collect(1:10) #[1,  4,  10]
paramsTempReg[:idxBG] =  push!(collect(1:2:(L-1)),L) # [1, 2, 5, 9,  10] #collect(1:10) #[1, 10]
cTempReg = reconstruction(fSF, f; paramsTempReg...).data.data


#ImageView.imshow(abs.(cat(cBGEst, cTempReg, cStaticCorr, dims=3)[1,:,:,1,:]),CLim(0,maximum(cTempReg)), name="TempReg")
ImageView.imshow(abs.(cat(cBGEst, cTRTikhonov, cTempReg,  dims=3)[1,:,:,1,:]),CLim(0,maximum(cBGEst)), name="TempReg")

#ImageView.imshow(abs.(cTempReg[1,:,:,1,:]),CLim(0,maximum(cTempReg)), name="TempReg")

# ###########################
# ## Visualize the results ##
# ###########################

#ImageView.imshow(abs.(cStaticCorr[1,:,:,1,:]),CLim(0,maximum(cStaticCorr)), name="Static")
#ImageView.imshow(abs.(cBGEst[1,:,:,1,:]),CLim(0,maximum(cBGEst)), name="Joint")


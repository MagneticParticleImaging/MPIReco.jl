using Pkg

# Install required packages
for P in ["HTTP","TensorDecompositions", "PyPlot"]
  !haskey(Pkg.installed(), P) && Pkg.add(P)
end

# Download data
include("downloadData.jl")

using MPIReco, DelimitedFiles, PyPlot
include("./utils/subsampling.jl")
include("./utils/FRTruncation.jl")

# load system matrix
datadirSM = "./data/"
bSF = MPIFile(datadirSM*"SM.mdf")
# backup for local testing
# bSF = MPIFile("/opt/mpidata/20180228_180048_OpenMPIData_1_1/9")
nx,ny,nz = calibSize(bSF);
f = [2145,8629,60903]
S = ComplexF64.( getSystemMatrix(bSF,f,bgCorrection=true, tfCorrection=false) )

# sampling locations (linear index)
samplingIdx = vec( readdlm("./data/pdSamplingIdx.txt", '\t', Int64) )

# get undersampled data
y = sfMeas(S, samplingIdx)

# CS-recovery
par_cs = Dict{Symbol,Any}()
par_cs[:shape] = (nx,ny,nz)   # size of the spatial measurement grid
par_cs[:λ_l1] = 0.02          # l1-regularization parameter
par_cs[:ρ_l1] = 2.0	      # parameter ρ associated with the l1-regularization
par_cs[:iterationsInner] = 50 # max. number of inner Split Bregman iterations
par_cs[:iterations] = 10      # number of outer Split Bregman iterations
par_cs[:relTol] = 1.e-2       # stopping criterion for the inner Split Bregman loop
par_cs[:absTol] = 1.e-4       # stopping criterion for the inner Split Bregman loop

S_cs = smRecovery(y, samplingIdx, par_cs)

# CSLR-recovery
par_cslr = Dict{Symbol,Any}()
par_cslr[:shape] = (nx,ny,nz)
par_cslr[:λ_l1] = 0.02
par_cslr[:ρ_l1] = 2.0
par_cslr[:λ_lr] = 0.02	        # TNN regularization parameter (for CSLR and CSFR)
par_cslr[:ρ_lr] = 0.2           # parameter ρ associated with the low rank regularization
par_cslr[:iterationsInner] = 50
par_cslr[:iterations] = 10
par_cslr[:relTol] = 1.e-2
par_cslr[:absTol] = 1.e-4

S_cslr = smRecovery(y, samplingIdx, par_cslr)

# reshape into tensor format
S = reshape(S,nx,ny,nz,:)
S_cs = reshape(S_cs,nx,ny,nz,:)
S_cslr = reshape(S_cslr,nx,ny,nz,:)

# CSFR
r = (6,6,6) # maximum HOSVD rank
S_csfr = FRTruncation3d(S_cslr,r)

# plot frequency components
figure()
subplot(4,3,1)
title("k=$(f[1])")
ylabel("measured")
imshow(abs.(S[:,:,17,1]),cmap="gray"); xticks([]); yticks([])
subplot(4,3,2)
title("k=$(f[2])")
imshow(abs.(S[:,:,17,2]),cmap="gray"); xticks([]); yticks([])
subplot(4,3,3)
title("k=$(f[3])")
imshow(abs.(S[:,:,17,3]),cmap="gray"); xticks([]); yticks([])
subplot(4,3,4)
ylabel("CS")
imshow(abs.(S_cs[:,:,17,1]),cmap="gray"); xticks([]); yticks([])
subplot(4,3,5)
imshow(abs.(S_cs[:,:,17,2]),cmap="gray"); xticks([]); yticks([])
subplot(4,3,6)
imshow(abs.(S_cs[:,:,17,3]),cmap="gray"); xticks([]); yticks([])
subplot(4,3,7)
ylabel("CSLR")
imshow(abs.(S_cslr[:,:,17,1]),cmap="gray"); xticks([]); yticks([])
subplot(4,3,8)
imshow(abs.(S_cslr[:,:,17,2]),cmap="gray"); xticks([]); yticks([])
subplot(4,3,9)
imshow(abs.(S_cslr[:,:,17,3]),cmap="gray"); xticks([]); yticks([])
subplot(4,3,10)
ylabel("CSFR")
imshow(abs.(S_csfr[:,:,17,1]),cmap="gray"); xticks([]); yticks([])
subplot(4,3,11)
imshow(abs.(S_csfr[:,:,17,2]),cmap="gray"); xticks([]); yticks([])
subplot(4,3,12)
imshow(abs.(S_csfr[:,:,17,3]),cmap="gray"); xticks([]); yticks([])

# # GPU Acceleration
include("../../download.jl") #hide
# MPIReco supports generic GPU acceleration. This means that the user can use any GPU array type that supports the GPUArrays interface. This includes [CUDA.jl](https://github.com/JuliaGPU/CUDA.jl), [AMDGPU.jl](https://github.com/JuliaGPU/AMDGPU.jl), and [Metal.jl](https://github.com/JuliaGPU/Metal.jl).
# To perform a reconstruction on the GPU, one has to load a GPU backend package such as CUDA and specify the GPU array type:
# ```julia
# using CUDA
# gpu = CuArray
# ```
gpu = Array; #hide # comment out when using proper GPUs

# Afterwards one can use the normal reconstruction interface and specify the `arrayType` parameter.
# However, the default solver is Kaczmarz, which has poor GPU performance. One should instead use CGNR or other solver.
# Since those solvers usually converge slower, one also has to increase the number of iterations
using MPIReco #hide
bSF = MPIFile(joinpath(datadir, "calibrations", "12.mdf")) #hide
b = MPIFile(joinpath(datadir, "measurements", "20211226_203916_MultiPatch", "1.mdf")) #hide
params = Dict{Symbol, Any}()
params[:SNRThresh] = 5
params[:frames] = 1:1
params[:minFreq] = 80e3
params[:recChannels] = 1:2
params[:spectralLeakageCorrection] = true
params[:sf] = bSF
params[:weightingParams] = WhiteningWeightingParameters(whiteningMeas = bSF)
params[:reg] = [L2Regularization(0.1f0)]

# Specifying the GPU array type to use, allows reconstruction algorithm to move the system matrix and data to the GPU after frequency filtering:
c = reconstruct("SinglePatch", b; params..., arrayType = gpu, iterations = 100, solver = CGNR);

# The resulting image is moved to the CPU at the end of the reconstruction and is accessible like any other reconstruction:
using CairoMakie
fig = heatmap(c[1, :, :, 1, 1].data.data)
hidedecorations!(fig.axis)
fig
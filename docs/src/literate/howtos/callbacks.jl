# # Callbacks
include("../../download.jl") #hide
using MPIReco #hide
bSF = MPIFile(joinpath(datadir, "calibrations", "12.mdf")) #hide
b = MPIFile(joinpath(datadir, "measurements", "20211226_203916_MultiPatch", "1.mdf")); #hide
# RegularizedLeastSquares provides a callback mechanism that allows you to access and monitor the state of the optimization process.
# These callbacks are invoked at the start of each iteration and have the form f(solver, iteration). Via the solver, you can access any internal property of the solver state.

# MPIReco exposes this interface for compatible algorithms. You can use do-syntax to pass a callback of the form f(solver, frame, iteration) to the reconstruction:
parameters = Dict{Symbol, Any}() #hide
parameters[:SNRThresh] = 5 #hide
parameters[:sf] = bSF #hide
parameters[:frames] = 1:acqNumFrames(b) #hide
parameters[:minFreq] = 80e3 #hide
parameters[:recChannels] = 1:rxNumChannels(b) #hide
parameters[:iterations] = 3 #hide
parameters[:spectralLeakageCorrection] = true; #hide
c = reconstruct("SinglePatch", b; parameters...) do solver, frame, iteration
  tmp = sum(solversolution(solver))
  @info "Sum of concentration at frame $frame and iteration $iteration is $tmp"
end

# MPIReco also provides built-in callbacks that can store the solution in each iteration or compare the current solution with a reference:
cb = StoreSolutionPerFrameCallback()
reconstruct(cb, "SinglePatch", b; parameters...);
cb.solutions[1]

# You can also combine multiple callbacks:
cb1 = StoreSolutionPerFrameCallback()
ref = reshape(c[1, :, :, :, 1].data.data, :, 1)
cb2 = CompareSolutionPerFrameCallback(ref)
reconstruct("SinglePatch", b; parameters...) do solver, frame, iteration
  cb1(solver, frame, iteration)
  cb2(solver, frame, iteration)
end
cb2.results[1]

# Note that callbacks are used directly in the solver and thus reflect its value domain. 
# This means that in sparse reconstruction, the current solution approximation is not in the image domain. Likewise, the solution might still be a complex number.

# Both variants shown above get passed as a callbacks keyword argument to the algorithm.
# When using the callback directly like this, the interface changes from f(solver, frame, iteration) to f(solver, iteration).
# And you have to manually differentiate between frames by watching for the start of a new reconstruction:
frame = 0
function my_callback(solver, iteration)
  if iteration == 0
    global frame += 1
  end
  @info "Frame $frame, Iteration $iteration"
end
reconstruct("SinglePatch", b; parameters..., callbacks = my_callback);
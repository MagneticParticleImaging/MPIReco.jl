export CompareSolutionPerFrameCallback, StoreSolutionPerFrameCallback
"""
    CompareSolutionPerFrameCallback(ref, cmp = nrmsd)

Callback that compares the solver's current `solution` with the given reference via `cmp(ref, solution)` for each iteration and frame. 
Results are stored in the `results` field as a vector of comparison results per frame.
"""
mutable struct CompareSolutionPerFrameCallback{T, F}
  ref::Matrix{T}
  cmp::F
  results::Vector{Vector{Float64}}
  frame::Int64
  CompareSolutionPerFrameCallback(ref::Matrix{T}, cmp::F = nrmsd) where {T, F} = new{T, F}(ref, cmp, [Float64[]], 1)
end
function (cb::CompareSolutionPerFrameCallback)(solver, frame, _)
  if frame != cb.frame
    push!(cb.results, Float64[])
    cb.frame = frame
  end
  c = solversolution(solver)
  push!(cb.results[frame], cb.cmp(cb.ref[:, frame], c))
end

"""
    StoreSolutionPerFrameCallback(T = Float32)

Callback that stores the solver's `solution` at each iteration, grouped by frame.
Results are stored in the `solutions` field. Aggregated solutions are a vector (frames) of vectors (iterations) of Vector{T} (solution data).
"""
mutable struct StoreSolutionPerFrameCallback{T}
  solutions::Vector{Vector{Vector{T}}}
  frame::Int64
  StoreSolutionPerFrameCallback(T::Type = Float32) = new{T}(Vector{Vector{Vector{T}}}(), 0)
end
function (cb::StoreSolutionPerFrameCallback{T})(solver, frame, _) where T
  if frame != cb.frame
    push!(cb.solutions, Vector{Vector{T}}())
    cb.frame = frame
  end
  c = deepcopy(solversolution(solver))
  push!(cb.solutions[frame], c)
end
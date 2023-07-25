export ThreadPinnedAlgorithm, ThreadPinnedAlgorithmParameter
Base.@kwdef struct ThreadPinnedAlgorithmParameter <: AbstractReconstructionAlgorithmParameter
  threadID::Int64
  algo::AbstractReconstructionAlgorithm
end

mutable struct ThreadPinnedAlgorithm <: AbstractReconstructionAlgorithm
  params::ThreadPinnedAlgorithmParameter
  recoTask::Union{Nothing,Task}
  taskLock::ReentrantLock
  inputChannel::Channel{Any}
  outputChannel::Channel{Any}
end

ThreadPinnedAlgorithm(params::ThreadPinnedAlgorithmParameter) = ThreadPinnedAlgorithm(params, nothing, ReentrantLock(), Channel{Any}(Inf), Channel{Any}(Inf))

take!(algo::ThreadPinnedAlgorithm) = take!(algo.outputChannel)
function put!(algo::ThreadPinnedAlgorithm, u)
  put!(algo.inputChannel, u)
  lock(algo.taskLock)
  try
    if isnothing(algo.recoTask) || istaskdone(algo.recoTask)
      algo.recoTask = @tspawnat algo.params.threadID pinnedRecoTask(algo)
    end
  finally
    unlock(algo.taskLock)
  end
end
function pinnedRecoTask(algo::ThreadPinnedAlgorithm)
  while isready(algo.inputChannel)
    result = nothing
    try
      put!(algo.params.algo, take!(algo.inputChannel))
      result = take!(algo.params.algo)
    catch e
      result = e
    end
    put!(algo.outputChannel, result)
  end
end
# TODO general async task, has to preserve order (cant just spawn task for each put)
# TODO Timeout task with timeout options for put and take
# TODO maybe can be cancelled?
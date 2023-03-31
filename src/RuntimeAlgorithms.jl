mutable struct ThreadPinnedReconstruction <: AbstractMPIReconstructionAlgorithm
  threadID::Int64
  recoTask::Union{Nothing,Task}
  taskLock::ReentrantLock
  algo::AbstractMPIReconstructionAlgorithm
  inputChannel::Channel{Any}
  outputChannel::Channel{Any}
end

take!(algo::ThreadPinnedReconstruction) = take!(algo.outputChannel)
function put!(algo::ThreadPinnedReconstruction, u)
  put!(algo.inputChannel, u)
  lock(algo.taskLock)
  try
    if isnothing(algo.recoTask) || istaskdone(algo.recoTask)
      algo.recoTask = @tspawnat algo.threadID pinnedRecoTask(algo)
    end
  finally
    unlock(algo.taskLock)
  end
end
function pinnedRecoTask(algo::ThreadPinnedReconstruction)
  while isready(algo.inputChannel)
    # TODO error handling
    put!(algo.algo, take!(algo.inputChannel))
    reco = take!(algo.algo)
    put!(algo.outputChannel, reco)
  end
end

outputType(algo::ThreadPinnedReconstruction) = IntermediateOutput()

# TODO general async task, has to preserve order (cant just spawn task for each put)
# TODO Timeout task with timeout options for put and take
# TODO maybe can be cancelled?
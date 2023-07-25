mutable struct Progress
  const lock::ReentrantLock
  descr::String
  counter::Int64
  const total::Int64

  Progress(descr::String, start::Int64, total::Int64) = new(ReentrantLock(), descr, start, total)
end


# What is total in an algorithm with nested algorithms
# How to handle multiple threads

function step!(prog::Progress, size::Int64=1; descr = description(prog))
  access!(prog) do
    prog.descr = descr
    prog.counter = min(prog.total, prog.counter + size)
    return prog.counter
  end
end

function state(prog::Progress)
  access!(prog) do
    return (prog.descr, prog.counter, prog.total) 
  end
end

function access!(f::Function, prog::Progress)
  lock(prog.lock)
  try
    f()
  finally
    unlock(prog.lock)
  end
end
using Requests

if !isdir("dataMP01")
  streamSM = get("http://media.tuhh.de/ibi/mpireco/data.zip")
  save(streamSM, "data.zip")
  run(`unzip data.zip`)
end

include("MultiPatch.jl")
include("Regular.jl")
include("Sparse.jl")

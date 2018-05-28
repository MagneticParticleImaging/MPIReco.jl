using Requests

if !isdir("dataMP01")
  streamSM = get("http://media.tuhh.de/ibi/mpireco/data.zip")
  save(streamSM, "data.zip")
  run(`unzip data.zip`)
end

filenameSM = "systemMatrix.mdf"
filenameMeas = "measurement.mdf"

if !isfile(filenameSM)
  streamSM = get("http://media.tuhh.de/ibi/mdfv2/systemMatrix_V2.mdf")
  save(streamSM, filenameSM)
end
if !isfile(filenameMeas)
  streamMeas = get("http://media.tuhh.de/ibi/mdfv2/measurement_V2.mdf")
  save(streamMeas, filenameMeas)
end

include("MultiPatch.jl")
include("MultiGradient.jl")
include("Regular.jl")
include("Sparse.jl")

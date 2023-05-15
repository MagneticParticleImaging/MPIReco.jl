module RecoUtils

using ThreadPools

import Base: put!, take!

include("AlgorithmInterface.jl")
include("StructTransforms.jl")
include("MiscAlgorithms/MiscAlgorithms.jl")

end # module

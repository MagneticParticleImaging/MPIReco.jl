include("Splines.jl")
include("Tikhonov.jl")

function reconstructionTempReg(params...; kargs...)
  if haskey(kargs, :λ) && length(kargs[:λ]) > 1 &&
     haskey(kargs, :β) && length(kargs[:β]) > 1
    reconstructionTempRegTikhonov(params...; kargs...)
  else
    reconstructionTempRegSplines(params...; kargs...)
  end
end


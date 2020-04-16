using TensorDecompositions

function FRTruncation3d(S::Array{T,4},r::NTuple{3,Int64}) where T
  S_lr = zeros(T,size(S));
  for k=1:size(S,4)
    S_lr[:,:,:,k] .= truncateHOSVD(S[:,:,:,k],r)
  end

  return S_lr
end

function truncateHOSVD(x::Array{T,3},r::NTuple{3,Int64}) where T
  guvw_r = hosvd(real.(x),r)
  guvw_i = hosvd(imag.(x),r)

  return compose(guvw_r)+1im*compose(guvw_i)
end


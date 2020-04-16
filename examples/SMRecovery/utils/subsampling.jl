function sfMeas(sf::Matrix{T}, idx::Vector{Int64}) where T<:Complex

  meas = zeros(ComplexF64, length(idx), size(sf,2))
  for k=1:size(sf,2)
    meas[:,k] .= sf[:,k][idx]
  end

  return meas
end


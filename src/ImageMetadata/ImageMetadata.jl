Base.@kwdef struct ImageAxisParameter <: AbstractMPIRecoParameters
  pixspacing::Vector{typeof(1.0u"mm")}
  offset::Vector{typeof(1.0u"mm")}
  dt::typeof(1.0u"s")
end

function RecoUtils.process(t::Type{AbstractMPIReconstructionAlgorithm}, data::Array, params::ImageAxisParameter)
  return makeAxisArray(data, params.pixspacing, params.offset, params.dt)
end

Base.@kwdef struct ImageMetadataSystemMatrixParameter <: AbstractMPIRecoParameters
  meas::MPIFile
  sm::Union{MPIFile, Vector{MPIFile}}
  grid::RegularGridPositions
  axis::ImageAxisParameter
end

function RecoUtils.process(t::Type{AbstractMPIReconstructionAlgorithm}, data::Array, params::ImageMetadataSystemMatrixParameter)
  numcolors = 1
  if isa(params.sm,AbstractVector) || isa(params.sm,MultiContrastFile)
    numcolors = length(params.sm)
  end
  
  shp = shape(params.grid)
  imArray = Array{Float32}(undef, numcolors, shp..., size(data)[end])
  im = process(t, imArray, params.axis)

  temp = reshape(data,reduce(*, shp),numcolors,:)
  temp = permutedims(temp, [2,1,3])
  im[:] = temp[:]

  return ImageMeta(im, generateHeaderDict(params.sm, params.meas))
end
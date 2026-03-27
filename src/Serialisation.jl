function StructUtils.lower(::MPIRecoStyle, file::MPIFile)
    return filepath(file)
end
function StructUtils.lift(::MPIRecoStyle, ::Type{T}, source::String) where {T<:MPIFile}
    return MPIFile(source), source
end
function StructUtils.lower(::MPIRecoStyle, file::MultiMPIFile)
    return filepath.(file)
end
function StructUtils.lift(::MPIRecoStyle, ::Type{MultiMPIFile}, source::Vector)
    files = MPIFile.(source)
    return MultiMPIFile(files), source
end
function StructUtils.lower(::MPIRecoStyle, file::MultiContrastFile)
    return filepath.(file)
end
function StructUtils.lift(::MPIRecoStyle, ::Type{MultiContrastFile}, source::Vector)
    files = MPIFile.(source)
    return MultiContrastFile(files), source
end


function StructUtils.lower(::MPIRecoStyle, range::UnitRange)
    return Dict{String, Any}(
        "start" => range.start,
        "stop" => range.stop
    )
end
function StructUtils.lift(::MPIRecoStyle, ::Type{T}, dict::Dict) where {T<:UnitRange}
    start = dict["start"]
    stop = dict["stop"]
    return UnitRange(start, stop), dict
end

function StructUtils.lower(::MPIRecoStyle, reg::T) where {T<:AbstractRegularization}
    dict = Dict{String, Any}(
        MODULE_TAG => string(parentmodule(T)),
        TYPE_TAG => string(nameof(T)),
        "λ" => reg.λ
    )
    
    # Add all other fields except λ
    for field in fieldnames(T)
        if field != :λ
            value = getfield(reg, field)
            # Skip empty sparse trafos or other empty collections
            if !(value isa AbstractArray && isempty(value))
                dict[string(field)] = StructUtils.lower(FIELD_STYLE[], value)
            end
        end
    end
    
    return dict
end

function StructUtils.lift(::MPIRecoStyle, ::Type{T}, dict::Dict{String, Any}) where {T<:AbstractRegularization}
    λ = dict["λ"]
    
    # Collect kwargs from all keys except λ and metadata
    filteredKeys = filter(x -> !(isequal("λ", x) || startswith(x, "_")), keys(dict))
    kwargs = Dict(Symbol(x) => dict[x] for x in filteredKeys)
    
    # Remove empty sparseTrafo if present
    if haskey(kwargs, :sparseTrafo) && isempty(kwargs[:sparseTrafo])
        pop!(kwargs, :sparseTrafo)
    end
    
    return T(λ; kwargs...), dict
end
function StructUtils.lift(style::MPIRecoStyle, ::Type{Vector{T}}, source::Vector) where {T<:AbstractRegularization}
    result = Any[]
    
    for val in source
        if val isa Dict && haskey(val, MODULE_TAG) && haskey(val, TYPE_TAG)
            # Has type metadata - get concrete type
            modDict = MODULE_DICT[]
            if isnothing(modDict)
                error("Cannot lift regularization without MODULE_DICT context")
            end
            
            mod_name = val[MODULE_TAG]
            type_name = val[TYPE_TAG]
            concrete_type = modDict[mod_name, type_name]
            
            if !isnothing(concrete_type)
                lifted, _ = StructUtils.lift(style, concrete_type, val)
                push!(result, lifted)
            else
                error("Cannot find type $type_name in module $mod_name")
            end
        else
            error("Regularization entry missing type metadata: $val")
        end
    end
    
    # Narrow vector to concrete type if possible
    return identity.(result), source
end

function StructUtils.lower(::MPIRecoStyle, norm::T) where {T<:AbstractRegularizationNormalization}
    return Dict{String, Any}(
        MODULE_TAG => string(parentmodule(T)),
        TYPE_TAG => string(nameof(T))
    )
end

function StructUtils.lift(::MPIRecoStyle, ::Type{T}, dict::Dict) where {T<:AbstractRegularizationNormalization}
    return T(), dict
end
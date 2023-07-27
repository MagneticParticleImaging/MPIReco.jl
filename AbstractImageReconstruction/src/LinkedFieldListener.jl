export LinkedFieldListener
struct LinkedFieldListener{T<:AbstractReconstructionAlgorithmParameter} <: SerializableListener
  plan::RecoPlan{T}
  field::Symbol
  fn::Function
end

function propertyupdate!(listener::LinkedFieldListener, origin::RecoPlan, name::Symbol, old::X, new::Y) where {X, Y}
  if hasmethod(listener.fn, Tuple{X, Y}) && !ispropertyset(listener.plan, listener.field)
    setvalue!(listener.plan, listener.field, listener.fn(old, new))
  end
end

function addDictValue!(dict, value::LinkedFieldListener)
  dict["plan"] = toDictValue(parentfields(value.plan))
  dict["field"] = toDictValue(value.field)
  dict["fn"] = toDict(value.fn)
end

function loadListener(::Type{LinkedFieldListener}, root::RecoPlan, dict, modDict)
  plan = root
  for param in dict["plan"][1:end]
    plan = plan[Symbol(param)]
  end
  field = Symbol(dict["field"])
  fn = tomlType(dict["fn"], modDict)
  return LinkedFieldListener(plan, field, fn)
end
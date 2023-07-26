export LinkedFieldListener
struct LinkedFieldListener{T<:AbstractReconstructionAlgorithmParameter} <: AbstractPlanListener
  plan::RecoPlan{T}
  field::Symbol
  fn::Function
end

function propertyupdate!(listener::LinkedFieldListener, origin::RecoPlan, name::Symbol, old::X, new::Y) where {X, Y}
  if hasmethod(listener.fn, Tuple{X, Y}) && !ispropertyset(listener.plan, listener.field)
    setvalue!(listener.plan, listener.field, listener.fn(old, new))
  end
end
"""
    TraitRelationship

The relationship between a trait and its environment, represented as a Matrix
of Floats.

"""
mutable struct TraitRelationship
  traitfun::Function
end

function GaussTemp(temp::Float64, opttemp::Float64, var::Float64)
  1/sqrt(2 * π * var^2) * exp(-abs(temp - opttemp)^2/(2 * var^2))
end

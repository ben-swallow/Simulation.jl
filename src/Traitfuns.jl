function traitfun(eco::Ecosystem, pos::Int64, spp::Int64)
    hab = eco.abenv.habitat
    trts = eco.spplist.traits
    rel = eco.relationship
  _traitfun(hab, trts, rel, pos, spp, iscontinuous(hab))
end
#function _traitfun(hab::AbstractHabitat, trts::AbstractTraits,
#     rel::AbstractTraitRelationship, pos::Int64,
#     spp::Int64, cont::Array{Bool, 1})
#    traitnames = fieldnames(typeof(trts))
#    habnames = fieldnames(typeof(hab))
#    relnames = fieldnames(typeof(rel))
#    results = map(1:length(traitnames)) do tr
#        thishab = gethabitat(hab, habnames[tr])
#        thistrt = getpref(trts, traitnames[tr])
#        thisrel = getrelationship(rel, relnames[tr])
#        _traitfun(thishab, thistrt, thisrel, pos, spp, cont[tr])
#    end
#    combineTR(rel)(results)
#end
function _traitfun(hab::HabitatCollection2, trts::TraitCollection2, rel::R, pos::Int64, spp::Int64, cont::Array{Bool, 1}) where R <: AbstractTraitRelationship
    res1 = _traitfun(hab.h1, trts.t1, rel.tr1, pos, spp, iscontinuous(hab.h1))
    res2 = _traitfun(hab.h2, trts.t2, rel.tr2, pos, spp, iscontinuous(hab.h2))
    return combineTR(rel)(res1, res2)
end
function _traitfun(hab::ContinuousHab, trts::GaussTrait,
    rel::R, pos::Int64, spp::Int64, cont::Bool) where R <: AbstractTraitRelationship
        h = gethabitat(hab, pos)
        mean, var = getpref(trts, spp)
    return rel(h, mean, var)
end
function _traitfun(hab::ContinuousTimeHab, trts::GaussTrait,
    rel::R, pos::Int64, spp::Int64, cont::Bool) where R <: AbstractTraitRelationship
        h = gethabitat(hab, pos)
        mean, var = getpref(trts, spp)
    return rel(h, mean, var)
end
function _traitfun(hab::ContinuousTimeHab, trts::TempBin,
    rel::R, pos::Int64, spp::Int64, cont::Bool) where R <: AbstractTraitRelationship
        h = gethabitat(hab, pos)
        (a, b, c, d) = getpref(trts, spp)
    return rel(Trapezoid(a, b, c, d), h)
end
function _traitfun(hab::ContinuousTimeHab, trts::RainBin,
    rel::R, pos::Int64, spp::Int64, cont::Bool) where R <: AbstractTraitRelationship
        h = gethabitat(hab, pos)
        (a, b) = getpref(trts, spp)
    return rel(Uniform(a, b), h)
end
function _traitfun(hab::DiscreteHab, trts::DiscreteTrait,
    rel::R, pos::Int64, spp::Int64, cont::Bool) where R <: AbstractTraitRelationship
        currentniche = gethabitat(hab, pos)
        preference = getpref(trts, spp)
    return rel(currentniche, preference)
end

#function TraitFun(eco::Ecosystem, pos::Int64, spp::Int64, cont::Bool = true)
#  T = gethabitat(eco, pos)
#  T_opt, Var = getpref(env, eco, spp)
#  return gettraitfun(eco)(T, T_opt, Var)
#end

#function TraitFun(eco::Ecosystem, pos::Int64, spp::Int64, cont::Bool = false)
#  T = gethabitat(eco, pos)
#  T_opt, Var = getpref(env, eco, spp)
#  return gettraitfun(eco)(T, T_opt, Var)
#end
#function TraitFun(env::Type{DiscreteHab{Niches}}, eco::Ecosystem, pos::Int64, spp::Int64)
#  currentniche = gethabitat(eco, pos)
#  preference = getpref(env, eco, spp)
#  return gettraitfun(eco)(currentniche, preference)
#end

#function TraitFun(env::Type{ContinuousHab{None}}, eco::Ecosystem, pos::Int64, spp::Int64)
#  return 1.0
#end
function getpref(traits::GaussTrait, spp::Int64)
  return traits.mean[spp], traits.var[spp]
end
function getpref(traits::TempBin, spp::Int64)
  return traits.dist[spp, :]
end
function getpref(traits::RainBin, spp::Int64)
  return traits.dist[spp, :]
end
function getpref(traits::DiscreteTrait, spp::Int64)
  return traits.val[spp]
end

function getpref(traits::T, field::Symbol) where T <: AbstractTraits
  return getfield(traits, field)
end

function getrelationship(rel::R, field::Symbol) where R <: AbstractTraitRelationship
  return getfield(rel, field)
end

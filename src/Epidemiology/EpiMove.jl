"""
    EpiMovement{MO <: AbstractMovement} <: AbstractMovement

Movement can happen at several different levels, local gaussian processes, `home`, and longer distance commutes, `work`.
"""
mutable struct EpiMovement{MO1 <: AbstractMovement, MO2 <: AbstractMovement} <: AbstractMovement
    home::MO1
    work::Union{Missing, MO2}
end

function EpiMovement(homekernels::Vector{K}, home_to_work::DataFrame) where K <: AbstractKernel
    home = AlwaysMovement{K, NoBoundary}(homekernels, NoBoundary())
    work = Commuting(home_to_work)
    return EpiMovement{typeof(home), typeof(work)}(home, work)
end

function EpiMovement(homekernels::Vector{K}) where K <: AbstractKernel
    home = AlwaysMovement{K, NoBoundary}(homekernels, NoBoundary())
    work = missing
    return EpiMovement{typeof(home), typeof(home)}(home, work)
end

mutable struct Commuting <: AbstractMovement
    home_to_work::DataFrame
end

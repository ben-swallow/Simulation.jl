"""
    EpiLandscape

Disease class abundances housed in the landscape. These are represented in both 2 dimensions (for computational efficiency in simulations) and 3 dimensions (to represent disease classes, their abundances and position in the grid).

"""
mutable struct EpiLandscape
  matrix::Matrix{Int64}
  matrix_v::Matrix{Int64}
  grid::AxisArray{Int64}
  seed::Vector{MersenneTwister}
  function EpiLandscape(human_abun::Matrix{Int64}, virus_abun::Matrix{Int64}, ax)
    d = length.(ax)
    grid = AxisArray(reshape(human_abun, d), ax)
    return new(human_abun, virus_abun, grid, [MersenneTwister(rand(UInt)) for _ in 1:Threads.nthreads()])
  end
  function EpiLandscape(human_abun::Matrix{Int64}, virus_abun::Matrix{Int64}, ax, d2::Tuple, seed::Vector{MersenneTwister})
    d = length.(ax)
    grid = AxisArray(reshape(human_abun, d), ax)
    return new(human_abun, virus_abun, grid, seed)
  end
end

virus(x::EpiLandscape) = x.matrix_v
human(x::EpiLandscape) = x.matrix

import Base.copy
function copy(gl::EpiLandscape)
    return EpiLandscape(copy(gl.matrix), size(gl.grid))
end

function Base.isapprox(gl_1::EpiLandscape, gl_2::EpiLandscape; kwargs...)
    return isapprox(gl_1.matrix, gl_2.matrix; kwargs...)
end

"""
    emptyepilandscape(epienv::GridEpiEnv, epilist::EpiList)

Function to create an empty EpiLandscape given a GridEpiEnv and a
EpiList.
"""
function emptyepilandscape(epienv::GridEpiEnv, epilist::EpiList)
  mat_human = zeros(Int64, counttypes(epilist.human, true), countsubcommunities(epienv))
  mat_virus = zeros(Int64, counttypes(epilist.virus, true), countsubcommunities(epienv))
  ax = AxisArrays.axes(epienv.initial_population)
  return EpiLandscape(mat_human, mat_virus, ax)
end

using Pkg

"""
    env_bool(key)

Checks for an enviroment variable and fuzzy converts it to a bool
"""
env_bool(key, default=false) = haskey(ENV, key) ? lowercase(ENV[key]) ∉ ["0","","false", "no"] : default


println("Adding SimulationData")
Pkg.add(PackageSpec(url="https://github.com/ScottishCovidResponse/SimulationData.jl"))
if !env_bool("TRAVIS_CI_BUILD")
    Pkg.build("SimulationData")
end

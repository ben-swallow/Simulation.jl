using Simulation
using Test
using Distributions
using Unitful.DefaultSymbols
using Simulation.Units

include("TestCases.jl")
epi = TestEpiSystemHousehold()

# Test number of individuals in households matches overall population size
@test sum(epi.abundances.matrix, dims = 2)[:, 1] == sum(epi.households.infection_status, dims = 1)[1, :]

# Test number of individuals in households matches population size in each grid square
for i in 1:size(epi.abundances.matrix, 2)
    @test sum(epi.abundances.matrix[:, i]) == sum(epi.households.infection_status[epi.households.gridID .== i, :])
end

# Alter population and check that this breaks the test
epi.abundances.matrix[1, :] .+= 1
@test sum(epi.abundances.matrix, dims = 2)[:, 1] != sum(epi.households.infection_status, dims = 1)[1, :]

# Reinstantiate households and test this has corrected the infection_status matrix
instantiate_households!(epi)
@test sum(epi.abundances.matrix, dims = 2)[:, 1] == sum(epi.households.infection_status, dims = 1)[1, :]

# Simulate forwards and check the same thing
simulate!(epi, 1month, 1day)
@test sum(epi.abundances.matrix, dims = 2)[:, 1] == sum(epi.households.infection_status, dims = 1)[1, :]
for i in 1:size(epi.abundances.matrix, 2)
    @test sum(epi.abundances.matrix[:, i]) == sum(epi.households.infection_status[epi.households.gridID .== i, :])
end

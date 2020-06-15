using Simulation
using Unitful
using Unitful.DefaultSymbols
using Simulation.Units
using Simulation.ClimatePref
using StatsBase
using Distributions
using AxisArrays
using HTTP

do_plot = false

# Download and read in population sizes for Scotland
file, io = mktemp()
r = HTTP.request("GET", "https://raw.githubusercontent.com/ScottishCovidResponse/temporary_data/master/human/demographics/scotland/data/demographics.h5")
write(io, r.body)
close(io)
scotpop = parse_hdf5(file, grid = "10km", component = "grid10km/10year/persons")

# Read number of age categories
age_categories = size(scotpop, 3)

# Set up simple gridded environment
area = (AxisArrays.axes(scotpop, 1)[end] + AxisArrays.axes(scotpop, 1)[2] -
    2 * AxisArrays.axes(scotpop, 1)[1]) *
    (AxisArrays.axes(scotpop, 2)[end] + AxisArrays.axes(scotpop, 2)[2] -
    2 * AxisArrays.axes(scotpop, 2)[1]) * 1.0

# Set population to initially have no individuals
abun_h = (
    Susceptible = fill(0, age_categories),
    Exposed = fill(0, age_categories),
    Asymptomatic = fill(0, age_categories),
    Presymptomatic = fill(0, age_categories),
    Symptomatic = fill(0, age_categories),
    Hospitalised = fill(0, age_categories),
    Recovered = fill(0, age_categories),
    Dead = fill(0, age_categories)
)
disease_classes = (
    susceptible = ["Susceptible"],
    infectious = ["Asymptomatic", "Presymptomatic", "Symptomatic"]
)
abun_v = (Virus = 0,)

disease_axis = Axis{:class}(collect(keys(abun_h)))
compartment_axes = (AxisArrays.axes(scotpop)..., disease_axis)
# Populate susceptibles according to actual population spread
initial_pop = cat(
    scotpop,
    zeros((size(scotpop)..., length(disease_axis)-1)),
    dims=ndims(scotpop)+1
)
initial_pop = AxisArray(initial_pop, compartment_axes);
initial_pop = permutedims(initial_pop, [3, 4, 1, 2])

susc = @view initial_pop[class=:Susceptible]
exposed = @view initial_pop[class=:Exposed]
# Weight samples by number of susceptibles
samples = sample(CartesianIndices(susc), weights(1.0 .* vec(susc)), 100)
for samp in samples
    susc[samp] > 0 || continue
    # Add to exposed
    exposed[samp] += 1
    # Remove from susceptible
    susc[samp] -= 1
end

# Set simulation parameters
numclasses = length(abun_h)
numvirus = length(abun_v)
birth_rates = fill(0.0/day, numclasses, age_categories)
death_rates = fill(0.0/day, numclasses, age_categories)
birth_rates[:, 2:4] .= uconvert(day^-1, 1/20years)
death_rates[1:end-1, :] .= uconvert(day^-1, 1/100years)
virus_growth_asymp = virus_growth_presymp = virus_growth_symp = fill(0.1/day, age_categories)
virus_decay = 1.0/day
beta_force = fill(10.0/day, age_categories)
beta_env = fill(10.0/day, age_categories)
ageing = fill(0.0/day, age_categories - 1) # no ageing for now

# Prob of developing symptoms
p_s = fill(0.96, age_categories)
# Prob of hospitalisation
p_h = fill(0.2, age_categories)
# Case fatality ratio
cfr_home = cfr_hospital = fill(0.1, age_categories)
# Time exposed
T_lat = 3days
# Time asymptomatic
T_asym = 5days
# Time pre-symptomatic
T_presym = 1.5days
# Time symptomatic
T_sym = 5days
# Time in hospital
T_hosp = 5days
# Time to recovery if symptomatic
T_rec = 11days

param = SEI3HRDGrowth(birth_rates, death_rates, ageing,
                      virus_growth_asymp, virus_growth_presymp, virus_growth_symp, virus_decay,
                      beta_force, beta_env, p_s, p_h, cfr_home, cfr_hospital,
                      T_lat, T_asym, T_presym, T_sym, T_hosp, T_rec)
param = transition(param, age_categories)

epienv = simplehabitatAE(298.0K, area, NoControl(), initial_pop)

# Dispersal kernels for virus and disease classes
dispersal_dists = fill(1.0km, numclasses * age_categories)
cat_idx = reshape(1:(numclasses * age_categories), age_categories, numclasses)
dispersal_dists[vcat(cat_idx[:, 3:5]...)] .= 20.0km
kernel = GaussianKernel.(dispersal_dists, 1e-10)
movement = AlwaysMovement(kernel)

# Traits for match to environment (turned off currently through param choice, i.e. virus matches environment perfectly)
traits = GaussTrait(fill(298.0K, numvirus), fill(0.1K, numvirus))
epilist = EpiList(traits, abun_v, abun_h, disease_classes,
                  movement, param, age_categories)
rel = Gauss{eltype(epienv.habitat)}()

# Create epi system with all information
epi = EpiSystem(epilist, epienv, rel)

# Run simulation
times = 2months; interval = 1day; timestep = 1day
Nsteps = floor(Int, times/timestep) + 1
abuns = zeros(Int64, size(epi.abundances.matrix)..., Nsteps)
@time simulate_record!(abuns, epi, times, interval, timestep)

if do_plot
    using Plots
    # View summed SIR dynamics for whole area
    display(plot_epidynamics(epi, abuns, classes=keys(abun_h)))
    display(plot_epiheatmaps(epi, abuns, steps = [21]))
end

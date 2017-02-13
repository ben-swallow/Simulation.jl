#Load packages
using Diversity
using Diversity.ShortNames
using Simulation
using Plots
using PhyloTrees
using Distributions
#Pkg.update()
# Create a partition
part=MatrixLandscape(reshape([1, 2, 3, 4, 5, 6, 7, 8], (2, 2, 2)), Habitats([1 2; 3 4]))

# Create an ecosystem
eco=Ecosystem(part, Species(), StringTraits(["A", "B"]))

# Calculate ordinariness of ecosystem
getordinariness!(eco)

b = β(eco)

sb = subdiv(β(eco), 3)
m = metadiv(β(eco), 3)

# Create a simple habitat matrix
mat=ones(10, 10)
# Populate habitat with 10,000 individuals from 50 species
LS=populate(50, 10000, Habitats(mat))
# Create ecosystem from habitat, making every species distinct and
# having the same trait
eco=Ecosystem(LS,Species(), StringTraits(repmat(["A"],50)))
# Create grid of species richness
sr=SR(eco, 10000)

# Plot
heatmap(LS.abundances[1, :, :])
heatmap(sr)

# Investigate alpha Diversity
alpha_bar=subdiv(ᾱ(eco),0)
heatmap(alpha_bar)

randtree=jtree(17, Exponential(0.1))
Plots.plot(randtree)

coaltree=jcoal(14, 5)
Plots.plot(coaltree)

tree=jcoal(17, 100)
Plots.plot(tree)
switch_rate=0.5
traits=["A","B", "C"]
trait_tree=assign_trait(tree,switch_rate, traits)
get_traits(trait_tree)
get_times(trait_tree)

Plots.plot(tree,markershape=:circle,
markercolor= [:blue,false],
markerstrokecolor=[:black,false, :red])


using StatsBase
using RCall
# Set up Habitat
species=50; individuals=10000
mat=create_habitat((10,10), ["A","B"], [0.4,0.6])

# Set up tree
tree=jcoal(50, 100)
assign_traits!(tree, 0.5, ["A","B"])
sp_trt=get_traits(tree, true)
# Create ecosystem
pop=populate(species, individuals, Niches(mat), sp_trt)
eco=Ecosystem(pop,Species(), StringTraits(sp_trt))

# Check species richness before
before=SR(eco, 100000)
@rlibrary("fields")
@rput before
R"image.plot(before)"
# Update by one timestep
update!(eco, 0.9, 0.6,0.5, 0.1)
# Check species richness after
after=SR(eco, 100000)
@rput after
R"image.plot(after)"


using StatsBase
using RCall
# p= amount of fragmentation, A = expected proportion of habitat
mat=random_habitat((100,100), ["A","B"], 0.5, [0.5,0.5])
hab= Array{Int64}(100,100)
hab[mat.=="A"]=1
hab[mat.=="B"]=2
@rput hab
R"
image(hab, legend = F, axes=F);
"

species=10; individuals=50000
# Set up tree
tree=jcoal(species, 100)
assign_traits!(tree, 0.2, ["A","B"])
sp_trt=get_traits(tree, true)
budg= Array{Float64}(50,50)
fill!(budg, 100)
energy=repmat([1], species)
# Try with a skewed distribution
pop=populate(species, individuals, Niches(mat), sp_trt, Budget(budg),
  Multinomial(individuals, rand(Dirichlet(species,1))))
eco=Ecosystem(pop, Species(), StringTraits(sp_trt), RealEnergy(energy))
maximum(mapslices(sum,eco.partition.abundances,1))

@rlibrary("fields")
a = subdiv(ᾱ(eco), 2)
@rput a
R"par(mfrow=c(1,2));image.plot(a,col=heat.colors(100), breaks=seq(0,max(a), length.out=101));image(hab, legend = F)"

p = subdiv(ρ̄(eco), 1)
@rput p
R"par(mfrow=c(1,2));image.plot(p,col=heat.colors(100), breaks=seq(0,max(p), length.out=101));image(hab, legend = F)"

b = subdiv(β(eco), 1)
@rput b
R"par(mfrow=c(1,2));image.plot(b,col=heat.colors(100), breaks=seq(0,max(b), length.out=101));image(hab, legend = F)"

g = subdiv(γ(eco), 1)
@rput g
R"par(mfrow=c(1,2));image.plot(b,col=heat.colors(100), breaks=seq(0,max(g), length.out=101));image(hab, legend = F)"

# Set up initial conditions
birth = 0.4
death = 0.4
move = 0.5
timestep = 1
# Check species richness before
before=subdiv(ᾱ(eco), 2)
@rlibrary("fields")
@rlibrary("grDevices")
@rput before
R"par(mfrow=c(1,2));image.plot(before,col=rainbow(50)[1:20], breaks=seq(0,20,1));image(hab, legend = F)"
R"pdf(file='Before_steady.pdf', paper='a4r',onefile=T,width=11.69,height=6)"
R"par(mfrow=c(1,2));image.plot(before,col=rainbow(20), breaks=seq(0,20,1));image(hab, legend = F)"
R"dev.off()"
for i in 1:10

# Update by one timestep
update!(eco, birth, death, move, timestep)
# Check species richness after
after=subdiv(ᾱ(eco), 2)
@rput after
R"par(mfrow=c(1,2));image.plot(before,col=rainbow(50)[1:20], breaks=seq(0,20,1));image(hab, legend = F)"
maximum(mapslices(sum,eco.partition.abundances,1))
@rput i
#R"par(mfrow=c(1,2));image.plot(after,col=rainbow(20), breaks=seq(0,50,1));image(hab, legend = F)"
R"pdf(file=paste('After_steady', i, '.pdf', sep=''), paper='a4r',onefile=T,width=11.69,height=6)"
R"par(mfrow=c(1,2));print(image.plot(after,col=rainbow(20), breaks=seq(0,20,1)));image(hab, legend = F)"
R"dev.off()"
end
R"par(mfrow=c(1,2));image.plot(before,col=rainbow(20), breaks=seq(0,20,1));image(hab, legend = F)"

#Simplest case - everything has same energy
numSpecies=3
numTraits=2
numNiches=2
energy_vec = [2, 3, 5]
sppl = SpeciesList(numSpecies, numTraits, Multinomial(5000, numSpecies),
                   energy_vec)
abenv = MatrixAbioticEnv(numNiches, (1,1), 50000)
eco = Ecosystem(sppl, abenv)

birth = 0.4
death = 0.4
move = 0.0
timestep = 0.1

times=1000
abun = zeros(times+1, numSpecies); ener = zeros(times+1)
abun[1,:] = eco.partition.abundances[:, 1]
ener[1] = sum(eco.spplist.abun .* eco.spplist.energy.energy)

for i in 1:times
update!(eco, birth, death, move, 1.0, 0.0, timestep)
abun[i+1, :] = eco.partition.abundances[: , 1]
ener[i+1] = sum(eco.partition.abundances[: , 1] .* eco.spplist.energy.energy)
end
#54000
abun=abun[1:1000,:];ener=ener[1:1000]
@rput abun
@rput ener
R"
 ymax=max(abun)+2000; par(xaxs='i',yaxs='i', mfrow=c(1,2));
 plot(abun[,1], type='l', ylim=c(0, ymax), ylab='Abundance', xlab='Time');
 for (i in 1:ncol(abun)){
 lines(abun[,i], col=i)};
 legend('topright', legend=1:ncol(abun), col=1:ncol(abun), pch=20);
 plot(ener, type='l', ylim=c(0, max(ener)+8000), ylab='Energy', xlab='Time')
"


using RCall
## Run new system of set up
numSpecies=3
numTraits=2
numNiches=2
energy_vec = [2, 3, 5]
sppl = SpeciesList(numSpecies, numTraits, Multinomial(5000, numSpecies),
                   energy_vec)
abenv = MatrixAbioticEnv(numNiches, (1,1), 50000)
eco = Ecosystem(sppl, abenv)

birth = 0.4
death = 0.4
move = 0.0
timestep = 0.1

times=5000
abun = zeros(times+1, numSpecies); ener = zeros(times+1)
abun[1,:] = eco.partition.abundances[:, 1]
ener[1] = sum(eco.spplist.abun .* eco.spplist.energy.energy)

for i in 1:times
update!(eco, birth, death, move, 1.0, 0.0, timestep)
abun[i+1, :] = eco.partition.abundances[: , 1]
ener[i+1] = sum(eco.partition.abundances[: , 1] .* eco.spplist.energy.energy)
end
#54000

@rput abun
@rput ener
R"pdf('One pop with energy cap.pdf', paper='a4r', onefile=T,width=10,height=8);
 ymax=max(abun)+2000; par(xaxs='i',yaxs='i', mfrow=c(1,2));
 plot(abun[,1], type='l', ylim=c(0, ymax), ylab='Abundance', xlab='Time');
 for (i in 1:ncol(abun)){
 lines(abun[,i], col=i)};
 legend('topright', legend=1:ncol(abun), col=1:ncol(abun), pch=20);
 plot(ener, type='l', ylim=c(0, max(ener)), ylab='Energy', xlab='Time');
 dev.off()
"
## For presentation, start with simple case where everything uses up 1 unit of energy

 ## TWO SQUARES
 using RCall
 ## Run new system of set up
 numSpecies=3
 numTraits=2
 numNiches=2
 energy_vec = [2, 3, 5]
 sppl = SpeciesList(numSpecies, numTraits, Multinomial(5000, numSpecies),
                    energy_vec)
 abenv = MatrixAbioticEnv(numNiches, (2,1), 50000)

 eco = Ecosystem(sppl, abenv)

 birth = 0.4
 death = 0.4
 move = 0.0
 timestep = 0.1

 gridSize= 2

 times=5000
 abun = zeros(times+1, numSpecies, gridSize); ener = zeros(times+1, gridSize)
 abun[1,:, 1] = eco.partition.abundances[:, 1]
 abun[1,:, 2] = eco.partition.abundances[:, 2]
 ener[1, 1] = sum(eco.partition.abundances[:,1] .* eco.spplist.energy.energy)
 ener[1, 2] = sum(eco.partition.abundances[:,2] .* eco.spplist.energy.energy)

 for i in 1:times
 update!(eco, birth, death, move, 1.0, 0.0, timestep)
 abun[i+1, :, 1] = eco.partition.abundances[: , 1]
 abun[i+1, :, 2] = eco.partition.abundances[: , 2]
 ener[i+1, 1] = sum(eco.partition.abundances[: , 1] .* eco.spplist.energy.energy)
 ener[i+1, 2] = sum(eco.partition.abundances[: , 2] .* eco.spplist.energy.energy)
 end
 #54000
pop1=abun[:,:,1]
pop2=abun[:,:,2]
 @rput pop1
 @rput pop2
 @rput ener
 R"pdf('Two pop.pdf', paper='a4r', onefile=T,width=10,height=8);
 ymax1=max(pop1);ymax2=max(pop2); par(xaxs='i',yaxs='i', mfrow=c(2,2));
  plot(pop1[,1], type='l', ylim=c(0, ymax1), ylab='Population 1', xlab='Time');
  for (i in 1:ncol(pop1)){
  lines(pop1[,i], col=i)};
  legend('topright', legend=1:ncol(pop1), col=1:ncol(pop1), pch=20);
  plot(pop2[,1], type='l', ylim=c(0, ymax2), ylab='Population 2', xlab='Time');
  for (i in 1:ncol(pop2)){
  lines(pop2[,i], col=i)};
  legend('topright', legend=1:ncol(pop2), col=1:ncol(pop2), pch=20);
  plot(ener[,1], type='l', ylab='Energy 1', xlab='Time');
  plot(ener[,2], type='l', ylab='Energy 2', xlab='Time');
  dev.off()
"

  ## FOUR SQUARES
  using RCall
  ## Run new system of set up
  numSpecies=3
  numTraits=2
  numNiches=2
  energy_vec = [2, 3, 5]
  sppl = SpeciesList(numSpecies, numTraits, Multinomial(5000, numSpecies),
                     energy_vec)
  abenv = MatrixAbioticEnv(numNiches, (2,2), 50000)

  eco = Ecosystem(sppl, abenv)

  birth = 0.4
  death = 0.4
  move = 0.01
  timestep = 0.1

  gridSize= 4

  times=1000
  abun = zeros(times+1, numSpecies, gridSize)
  for g in 1:gridSize
  abun[1,:, g] = eco.partition.abundances[:, g]
  end

  for i in 1:times
    for j in 1:gridSize
    update!(eco, birth, death, move, 1.0, 0.0, timestep)
    abun[i+1, :, j] = eco.partition.abundances[: , j]
    end
  end
  #54000
 pop1=abun[:,:,1]
 pop2=abun[:,:,2]
 pop3=abun[:,:,3]
 pop4=abun[:,:,4]
  @rput pop1
  @rput pop2
  @rput pop3
  @rput pop4
  R"
   ymax1=max(pop1);ymax2=max(pop2);  ymax3=max(pop3);ymax4=max(pop4);
   par(xaxs='i',yaxs='i', mfrow=c(2,2));
   plot(pop1[,1], type='l', ylim=c(0, ymax1), ylab='Population 1', xlab='Time');
   for (i in 1:ncol(pop1)){
   lines(pop1[,i], col=i)};
   legend('topright', legend=1:ncol(pop1), col=1:ncol(pop1), pch=20);
   plot(pop2[,1], type='l', ylim=c(0, ymax2), ylab='Population 2', xlab='Time');
   for (i in 1:ncol(pop2)){
   lines(pop2[,i], col=i)};
   legend('topright', legend=1:ncol(pop2), col=1:ncol(pop2), pch=20);
   plot(pop3[,1], type='l', ylim=c(0, ymax3), ylab='Population 3', xlab='Time');
   for (i in 1:ncol(pop3)){
   lines(pop3[,i], col=i)};
   legend('topright', legend=1:ncol(pop3), col=1:ncol(pop3), pch=20);
   plot(pop4[,1], type='l', ylim=c(0, ymax4), ylab='Population 4', xlab='Time');
   for (i in 1:ncol(pop4)){
   lines(pop4[,i], col=i)};
   legend('topright', legend=1:ncol(pop4), col=1:ncol(pop4), pch=20);

"


#pdf('Four populations with migration & traits.pdf', paper='a4', onefile=T,width=8,height=10);
#   dev.off()

module Mating
include("Initialisation.jl")
using Main.Initialisation
using StatsBase: sample, Weights
using Random: shuffle

export local_mating, global_mating
function calculate_fitness(individual::Individual, location::Int64,selection_cofficient::Float64)::Float64
    #returns:
        #1 if the individual has no deliterious alleles
        #1 - selection_cofficient if the individual has one deliterious allele
        #1 - 2*selection_cofficient if the individual has two deliterious alleles
    deliterios_alleles = (location < number_of_patches/2) ? sum(individual.genome) : sum(broadcast(~,individual.genome)) #if the individual is in the first half of the patches, the first allele is the deliterious one, otherwise the second allele is the deliterious one
    return 1 - selection_cofficient*deliterios_alleles
end
#this mating function pools all gametes in one pool
function local_mating(patch::Patch, selection_cofficient::Float64)::Patch
    location = patch.location
    #select parents based on fitness
    male_fitnesses, female_fitnesses = ((male -> calculate_fitness(male,location, selection_cofficient)).(patch.males), (female -> calculate_fitness(female,location, selection_cofficient)).(patch.females))
    male_parents = sample(patch.males, Weights(male_fitnesses), individuals_per_patch)
    female_parents = sample(patch.females, Weights(female_fitnesses), individuals_per_patch)

    #create the offspring genomes
    male_gametes = (individual -> individual.genome[rand([1,2])]).(male_parents)
    female_gametes = (individual -> individual.genome[rand([1,2])]).(female_parents)
    offspring_genomes = zip(female_gametes, shuffle(male_gametes))|>collect
    female_genomes = offspring_genomes[1:Int(floor(individuals_per_patch/2))]
    male_genomes = offspring_genomes[Int(floor(individuals_per_patch/2)+1):end]

    #create the offspring, id is unique for each individual in the population
    female_offspring = [Individual(genome, false, patch.location*individuals_per_patch/2*1 + i) for (i, genome) in enumerate(female_genomes)]
    male_offspring = [Individual(genome, true, patch.location*individuals_per_patch/2*2 + i) for (i, genome) in enumerate(male_genomes)]

    #return the new patch
    return Patch(female_offspring,male_offspring, patch.location)
end

#Apply the local mating function to all patches in the population
function global_mating(population::Vector{Patch}, selection_cofficient::Float64)::Vector{Patch}
    local_mating_partial = ((patch) -> local_mating(patch, selection_cofficient))
    return local_mating_partial.(population)
end
end

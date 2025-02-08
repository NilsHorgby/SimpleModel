module Mating

using .Main.Initialisation
using StatsBase: sample, Weights
using Random: shuffle

export local_mating, global_mating

function calculate_fitness(individuals::BitMatrix, location::Int64,selection_cofficient::Float64, number_of_patches::Int64)::Vector{Float64}
    #input:
        #individuals: a BitMatrix of size (individuals_per_patch,2)
        #location: the location of the patch
        #selection_cofficient: the selection cofficient
        #number_of_patches: the number of patches
    #returns:
        #1 if the individual has no deliterious alleles
        #1 - selection_cofficient if the individual has one deliterious allele
        #1 - 2*selection_cofficient if the individual has two deliterious alleles
    deliterios_alleles = (location <= number_of_patches/2) ? sum(eachcol(individuals)) : sum(eachcol(.~individuals)) #if the individual is in the first half of the patches, the first allele is the deliterious one, otherwise the second allele is the deliterious one
    return 1 .- selection_cofficient .* deliterios_alleles
end

function random_allele(individuals::BitMatrix)::BitVector
    individuals_per_patch = size(individuals,1)
    indices = collect(1:individuals_per_patch) .* rand(1:2,individuals_per_patch)
    return individuals[indices]
end

#this mating function pools all gametes in one pool
function local_mating(patch::Patch, selection_cofficient::Float64, number_of_patches::Int64, individuals_per_patch::Int64)::Patch
    location = patch.location
    #select parents based on fitness
    male_fitnesses, female_fitnesses = calculate_fitness.([patch.males,patch.females],location, selection_cofficient, number_of_patches)
    #there are 2*N parents (in total), as we want N offspring
    male_parents = sample(eachrow(patch.males), Weights(male_fitnesses), individuals_per_patch) |> x -> reduce(vcat,transpose(x))
    female_parents = sample(eachrow(patch.females), Weights(female_fitnesses), individuals_per_patch) |> x -> reduce(vcat,transpose(x))

    #create the offspring genomes
    male_gametes, female_gametes = random_allele(male_parents), random_allele(female_parents)
    
    offspring_genomes = hcat(male_gametes, shuffle(female_gametes)) 
    females = offspring_genomes[1:Int(floor(individuals_per_patch/2)),:]
    males = offspring_genomes[Int(floor(individuals_per_patch/2)+1):end,:]

    #return the new patch
    return Patch(females, males, patch.location)
end

#Apply the local mating function to all patches in the population
function global_mating(population::Vector{Patch}, selection_cofficient::Float64, number_of_patches::Int64, individuals_per_patch::Int64)::Vector{Patch}
    return local_mating.(population, selection_cofficient,number_of_patches, individuals_per_patch)
end

global_mating(create_population(10,10), 0.1, 10, 10)


end


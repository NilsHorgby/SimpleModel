module Initialisation

export individuals_per_patch, number_of_patches, Individual, Patch, create_individual, create_subpopulation, create_population
#Constants
const individuals_per_patch = 100
const number_of_patches = 100


#An individual has a genome, a sex and an id
struct Individual
    genome:: Tuple{Bool, Bool}
    sex:: Bool #false = female, true = male
    id:: Int
end

#A patch has a location and a list of individuals
struct Patch
    females:: Vector{Individual}
    males:: Vector{Individual}
    location:: Int
end


#create the initial population
function create_individual(id,sex)
    genome = (rand(Bool), rand(Bool))
    return Individual(genome,sex,id)
end

#create the population of a single patch
function create_subpopulation(location)
    females = [create_individual(location*individuals_per_patch/2*1 + i,true) for i in 1:Int(floor(individuals_per_patch/2))]
    males = [create_individual(location*individuals_per_patch/2*2 + i,false) for i in 1:Int(floor(individuals_per_patch/2))]
    return Patch(females,males,location)
end

function create_population()
    patches = [create_subpopulation(i) for i in 1:number_of_patches]
    return patches
end
end
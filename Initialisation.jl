module Initialisation
using Random: bitrand
export Individual, Patch, create_individual, create_subpopulation, create_population


#A patch has a location and a list of individual-genomes
struct Patch
    females:: BitMatrix
    males:: BitMatrix
    location:: Int
end




#create the population of a single patch
function create_subpopulation(location::Int64, individuals_per_patch::Int64)
    females = bitrand((individuals_per_patch,2))
    males = bitrand((individuals_per_patch,2))
    return Patch(females,males,location)
end

function create_population(number_of_patches::Int64, individuals_per_patch::Int64)
    return create_subpopulation.(1:number_of_patches,individuals_per_patch)
end
end


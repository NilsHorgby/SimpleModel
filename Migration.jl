module Migration

using .Main.Initialisation
using Distributions: Normal

export migration

struct SexSpecificPatch
    individuals::Vector{Individual}
    sex::Bool
    location::Int
end

function put_within_boundaries(location::Int64,signed_distance::Int64, number_of_patches::Int64)
    new_patch::Int64 = location + signed_distance
    if new_patch <= 0
        return 1
    elseif new_patch > number_of_patches
        return number_of_patches
    else 
        return new_patch
    end
end

function discrete_distrubution(dispersal_distance::Float64)
    return Int(round(rand(Normal(0,dispersal_distance),1)[1]))
end

function migration(population::Vector{Patch},dispersal_distance::Float64, number_of_patches::Int64)::Vector{Patch}
    #population = deepcopy(population) the population will be modified in place, but this doesn't matter
    current_females = [SexSpecificPatch(patch.females,false, patch.location) for patch in population]
    current_males = [SexSpecificPatch(patch.males,true, patch.location) for patch in population]
    next_females::Vector{SexSpecificPatch} = [SexSpecificPatch([],false,population[i].location)
                                              for i=eachindex(population)]
    next_males::Vector{SexSpecificPatch} = [SexSpecificPatch([],true,population[i].location)
                                            for i=eachindex(population)]
    #migration one by one:
    
    function migrate_one_individual(population::Vector{SexSpecificPatch},next_population::Vector{SexSpecificPatch}, dispersal_distance::Float64)
        individual::Individual= population[1].individuals[1]
        signed_distance::Int64 = discrete_distrubution(dispersal_distance)
        deleteat!(population[1].individuals,1)
        
        #Boundary condition: if the individual tries to go outside the area, it doesn't move
        new_patch = put_within_boundaries(population[1].location,signed_distance, number_of_patches)
        append!(next_population[new_patch].individuals,[individual])  
        if length(population[1].individuals) == 0
            deleteat!(population,1)
        end
        return population,next_population
    end
    
    #we should be able to do tail recursion here
    function migrate(population::Vector{SexSpecificPatch},next_population::Vector{SexSpecificPatch}, dispersal_distance::Float64)
        while population != []
            population, next_population = migrate_one_individual(population,next_population,dispersal_distance)
        end
        return next_population
    end

    
    function combine_males_and_females(females::SexSpecificPatch,males::SexSpecificPatch)::Patch
        @assert females.location == males.location "The locations are not the same"
        return Patch(females.individuals,males.individuals, females.location)
    end
    population = combine_males_and_females.(migrate(current_females,next_females,dispersal_distance),migrate(current_males,next_males,dispersal_distance))
    return population
end

end
"""
Benchmark for migration function

BenchmarkTools.Trial: 10000 samples with 1 evaluation.
 Range (min … max):  227.133 μs …  17.273 ms  ┊ GC (min … max):  0.00% … 90.34%
 Time  (median):     266.618 μs               ┊ GC (median):     0.00%
 Time  (mean ± σ):   387.354 μs ± 827.794 μs  ┊ GC (mean ± σ):  25.67% ± 11.62%

  █▃▂▁                                                          ▁
  ████▆▅▁▁▁▁▃▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▄▅▅▆▅▆▆▇ █
  227 μs        Histogram: log(frequency) by time       5.91 ms <

 Memory estimate: 1.14 MiB, allocs estimate: 21494.
"""

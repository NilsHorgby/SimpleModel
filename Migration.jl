module Migration
include("Initialisation.jl")
using .Initialisation
#using .Main.Initialisation
using Distributions: Normal

export migration

function put_within_boundaries(new_patch::Int64, number_of_patches::Int64)
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


array = bitrand((10,2))

array[1,:] = [0,1]
array[1,:]

function migration(population::Vector{Patch},dispersal_distance::Float64)::Vector{Patch}
    number_of_patches::Int64 = length(population)
    number_individuals::Int64 = size(population[1].females)[1]
    max_number_individuals_per_patch = round(Int,(1 + dispersal_distance) * 4 * size(population[1].females)[1])
    current_females::Vector{BitMatrix} = [patch.females for patch in population]
    current_males::Vector{BitMatrix} = [patch.males for patch in population]
    next_females::Vector{BitMatrix} = [BitMatrix(undef,max_number_individuals_per_patch,2) for i=eachindex(population)]
    next_males::Vector{BitMatrix} = [BitMatrix(undef,max_number_individuals_per_patch,2) for i=eachindex(population)]
    
    
    function migrate(population::Vector{BitMatrix},next_population::Vector{BitMatrix}, dispersal_distance::Float64)
        #signed_distances::Matrix{Int64} =  round.(Int,rand(Normal(0,dispersal_distance),(number_of_patches, number_of_in_patches)))
        #next_locations::Matrix{Int64} = put_within_boundaries.(signed_distances .+ hcat([collect(1:number_of_patches) for i in 1:number_individuals]...),number_of_patches)
        number_of_inds_in_patches = zeros(Int,number_of_patches)
        @inbounds for location in 1:number_of_patches
            for individual_index in 1:number_individuals
                next_location::Int64 = put_within_boundaries(round(Int,rand(Normal(location,dispersal_distance))), number_of_patches)
                number_of_inds_in_patches[next_location] += 1
                next_index::Int64 = number_of_inds_in_patches[next_location]
                next_population[next_location][next_index,:] = population[location][individual_index,:]
            end
        end
        #((patch,number_of_inds_in_patch) -> patch[1:number_of_inds_in_patch,:]).(next_females,number_of_inds_in_patches )
        return [patch[1:number_of_inds_in_patch,:] for (patch,number_of_inds_in_patch) in zip(next_females,number_of_inds_in_patches)]
    end

    
    return Patch.(migrate(current_females,next_females,dispersal_distance),migrate(current_males,next_males,dispersal_distance),1:number_of_patches)
end

using BenchmarkTools
population = create_population(100,100)
population[1].males
migration(population,1.5)
max(length.((x -> x.males).(migration(population,1.5)))...)
migration(population,1.5)[1].males
end
"""
BenchmarkTools.Trial: 6262 samples with 1 evaluation.
Range (min … max):  541.119 μs …   7.071 ms  ┊ GC (min … max):  0.00% … 88.25%
Time  (median):     648.125 μs               ┊ GC (median):     0.00%
Time  (mean ± σ):   796.221 μs ± 750.201 μs  ┊ GC (mean ± σ):  18.90% ± 16.57%

▆█▆▃                                                        ▁ ▁
████▇▅▁▃▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▃▆▇▇▇▇███ █
541 μs        Histogram: log(frequency) by time       4.58 ms <

Memory estimate: 1.91 MiB, allocs estimate: 60613.
"""

module Migration
#include("Initialisation.jl")
#using .Initialisation
using .Main.Initialisation
using Distributions: Normal

export migration

#const max_range::Int64 = 5
function put_within_boundaries(new_patch::Int64, number_of_patches::Int64)
    if new_patch <= 0
        return 1
    elseif new_patch > number_of_patches
        return number_of_patches
    else 
        return new_patch
    end
end
#using SpecialFunctions: erf

function discrete_distrubution(dispersal_distance::Real, location::Int64)::Int64
    return round(Int,rand(Normal(location,dispersal_distance)))
end

#function discrete_distrubution_pdf(dispersal_distance:: Real)::Vector{Float64}
#    return [erf((i+0.5)/(dispersal_distance*sqrt(2))) - erf((i-0.5)/(dispersal_distance*sqrt(2))) for i in -max_range:max_range]
#end


function migrate(population::Vector{Population},next_population::Vector{Population}, dispersal_distance::Float64)::Vector{Vector{Tuple{Bool,Bool}}}
    number_of_patches::Int64 = length(population)
    @inbounds for location in 1:number_of_patches
        @inbounds for individual_index in 1:length(population[location]) #see if a while loop is faster
            next_location::Int64 = put_within_boundaries(discrete_distrubution(dispersal_distance,location), number_of_patches)
            push!(next_population[next_location],population[location][individual_index])
        end
    end
    return next_population
end 





function migration(population::Vector{Patch}, dispersal_distance::Float64)::Vector{Patch}
    number_of_patches::Int64 = length(population)
    #max_number_individuals_per_patch::Int64 = 250 #What the optimum here is a good question
    current_females::Vector{Population} = [patch.females for patch in population]
    current_males::Vector{Population} = [patch.males for patch in population]
    next_females::Vector{Population} = [Vector{Tuple{Bool,Bool}}() for i=eachindex(population)]
    next_males::Vector{Population} = [Vector{Tuple{Bool,Bool}}() for i=eachindex(population)]
    return Patch.(migrate(current_females,next_females,dispersal_distance),migrate(current_males,next_males,dispersal_distance),1:number_of_patches)
end


end

"""
BenchmarkTools.Trial: 10000 samples with 1 evaluation.
 Range (min … max):  122.934 μs …   6.500 ms  ┊ GC (min … max): 0.00% … 92.94%
 Time  (median):     139.813 μs               ┊ GC (median):    0.00%
 Time  (mean ± σ):   153.675 μs ± 142.818 μs  ┊ GC (mean ± σ):  4.96% ±  5.89%

          ▃██▆▅▂                                                 
  ▁▂▂▁▂▃▅████████▇▆▅▅▄▄▃▃▃▃▃▃▄▄▄▅▅▄▄▃▃▂▂▂▂▂▂▂▂▁▁▁▁▁▁▂▁▁▁▁▁▁▁▁▁▁ ▃
  123 μs           Histogram: frequency by time          199 μs <

 Memory estimate: 121.69 KiB, allocs estimate: 811.

"""
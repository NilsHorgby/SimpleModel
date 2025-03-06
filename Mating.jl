module Mating

using .Main.Initialisation
using StatsBase: sample, Weights, counts
using Random: shuffle!
using Distributions: Poisson, Binomial
using Statistics: mean

export local_mating, global_mating, poisson_mating, multinomial_mating




function poisson_gamete_production(individuals::Vector{Tuple{Bool,Bool}}, number_of_inds_in_patch::Int64, location::Int64, fittnesses_by_genotype_first_half::Tuple{Float64,Float64,Float64}, number_of_patches::Int64, carring_capacity::Int64,reproductive_rate::Float64)::BitVector
    #this function produces the gametes of the next generation
    #the "middle" of the habiat is at number_of_patches ÷ 2 + 0.5
    fittnesses_by_genotype::Tuple{Float64,Float64,Float64} = (location <= number_of_patches/2) ? fittnesses_by_genotype_first_half :
                                                                                                 reverse(fittnesses_by_genotype_first_half)
    individual_genotypes::Vector{Int64} = sum.(individuals) 
    n_AA::Int64 = sum(individual_genotypes .== 0) 
    n_Aa::Int64 = sum(individual_genotypes .== 1)  
    n_aa::Int64 = length(individual_genotypes) - n_AA - n_Aa

    mean_offspring_per_genotype::Tuple{Float64,Float64,Float64} = 2*exp(reproductive_rate*(1 - number_of_inds_in_patch/carring_capacity)) .* fittnesses_by_genotype
    mean_offspring::Float64 = sum(mean_offspring_per_genotype .* (n_AA,n_Aa,n_aa))
    relized_number_of_offspring::Int64 = rand(Poisson(mean_offspring))
    fittness_adjusted_ps::Tuple{Float64,Float64,Float64} = fittnesses_by_genotype .* (n_AA, n_Aa, n_aa)
    p::Float64 = (fittness_adjusted_ps[1] + fittness_adjusted_ps[2]*0.5) / sum(fittness_adjusted_ps)

    return rand(relized_number_of_offspring) .> p 
end 



function poisson_mating(patch::Patch, fittnesses_by_genotype_first_half::Tuple{Float64,Float64,Float64}, number_of_patches::Int64, carring_capacity::Int64,reproductive_rate::Float64)::Vector{Tuple{Bool,Bool}}
    number_of_inds_in_patch::Int64 = length(patch.males) + length(patch.females)
    male_gametes::BitVector  = poisson_gamete_production(patch.males, number_of_inds_in_patch, patch.location, fittnesses_by_genotype_first_half, number_of_patches, carring_capacity, reproductive_rate)
    female_gametes::BitVector = poisson_gamete_production(patch.females, number_of_inds_in_patch, patch.location, fittnesses_by_genotype_first_half, number_of_patches, carring_capacity, reproductive_rate)
    number_offspring::Int64 = min(length(male_gametes),length(female_gametes))
    return tuple.(male_gametes[1:number_offspring],female_gametes[1:number_offspring])
end

function multinomial_gamete_production(individuals::Vector{Tuple{Bool,Bool}}, location::Int64, fittnesses_by_genotype_first_half::Tuple{Float64,Float64,Float64}, number_of_patches::Int64, individuals_per_patch::Int64)::BitVector
    fittnesses_by_genotype::Tuple{Float64,Float64,Float64} = (location <= number_of_patches ÷ 2) ? fittnesses_by_genotype_first_half :
                                                                                                   reverse(fittnesses_by_genotype_first_half)
    n_AA::Int64 = sum(sum.(individuals) .== 0)
    n_aA::Int64 = sum(sum.(individuals) .== 1)
    n_aa::Int64 = individuals_per_patch - n_AA - n_aA

    ps::Tuple{Float64,Float64,Float64} = fittnesses_by_genotype .* (n_AA, n_aA, n_aa)
    p::Float64 = (ps[1] + 0.5*ps[2]) /sum(ps)
    return rand(2*individuals_per_patch) .> p
end

function multinomial_mating(patch::Patch, fittnesses_by_genotype_first_half::Tuple{Float64,Float64,Float64}, number_of_patches::Int64, individuals_per_patch::Int64)::Vector{Tuple{Bool,Bool}}
    male_gametes::BitVector = multinomial_gamete_production(patch.males, patch.location, fittnesses_by_genotype_first_half, number_of_patches, individuals_per_patch)
    female_gametes::BitVector = multinomial_gamete_production(patch.females, patch.location, fittnesses_by_genotype_first_half, number_of_patches, individuals_per_patch)
    return tuple.(male_gametes,female_gametes)
end

#this mating function pools all gametes in one pool
#the gamete_producing_function should be a partial function only having individuals as input
#ex: (ind -> binomial_gamete_production(ind, location, selection_cofficient, number_of_patches))
function local_mating(patch::Patch, mating_function::Function)::Patch
    #the gamete producing function must take a Vector{Tuple{Bool,Bool}} as input
    @assert hasmethod(mating_function, Tuple{Patch}) "The mating does not implement the correct method"

    #If there are excess gametes of one sex, randomly discard the excess
    offspring_genomes::Vector{Tuple{Bool,Bool}} = mating_function(patch)
    
    #Half of the offspring become males and half become females
    individuals_in_patch::Int64 = length(offspring_genomes)
    females::Vector{Tuple{Bool,Bool}} = offspring_genomes[1:floor(Int,individuals_in_patch/2)]
    males::Vector{Tuple{Bool,Bool}} = offspring_genomes[floor(Int,individuals_in_patch/2)+1:end]

    return Patch(females, males, patch.location)
end

#Apply the local mating function to all patches in the population
# selection_cofficient::Float64, individuals_per_patch::Int64, carring_capacity::Int64, reproductive_rate::Float64
function global_mating(population::Vector{Patch}, mating_function::Function)::Vector{Patch}
    return local_mating.(population, mating_function)
end

#example usage:
#if isinteractive() #only runs if the file is run in REPL
#selection_cofficient = 0.3
#number_of_patches = 100
#individuals_per_patch = 100
#gamete_producing_function = ((Inds, location)-> binomial_gamete_production(Inds, location, selection_cofficient, number_of_patches, individuals_per_patch))
#
#
#
#gamete_producing_function(population[1].males,1)
#population = create_population(100,100)
#global_mating(population,gamete_producing_function)
#end
end#module
"""
BenchmarkTools.Trial: 10000 samples with 1 evaluation.
 Range (min … max):  249.731 μs …   6.898 ms  ┊ GC (min … max):  0.00% … 95.43%
 Time  (median):     279.159 μs               ┊ GC (median):     0.00%
 Time  (mean ± σ):   336.582 μs ± 272.498 μs  ┊ GC (mean ± σ):  13.92% ± 14.64%

  ▇█▆▄▂   ▁                                                     ▂
  ██████▇████▇▅▄▄▁▄▃▃▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▄▃▄▅▆▇█▇████ █
  250 μs        Histogram: log(frequency) by time       1.71 ms <

 Memory estimate: 1.00 MiB, allocs estimate: 6203.
"""


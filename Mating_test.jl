using Distributions
using Random: bitrand
using StatsBase: shuffle!
cd("/home/gadus/programing-projects/julia-projects/SimpleModel")
include("Initialisation.jl")
include("Mating.jl")
include("Migration.jl")
include("Utils.jl")
using .Initialisation, .Utils
using Statistics: mean,median, std
using Distributions: Normal
using JLD2
using Plots
function poisson_gamete_production1(individuals::Vector{Tuple{Bool,Bool}}, number_of_inds_in_patch::Int64, location::Int64, fittnesses_by_genotype_first_half::Vector{Float64}, number_of_patches::Int64, carring_capacity::Int64,reproductive_rate::Float64)::BitVector
    #this function produces the gametes of the next generation
    #the "middle" of the habiat is at number_of_patches ÷ 2 + 0.5
    fittnesses_by_genotype::Vector{Float64} = (location <= number_of_patches ÷ 2) ? fittnesses_by_genotype_first_half :
                                                                                    reverse(fittnesses_by_genotype_first_half)
    individual_genotypes::Vector{Int64} = sum.(individuals) 
    n_AA::Int64 = sum(individual_genotypes .== 0) 
    n_Aa::Int64 = sum(individual_genotypes .== 1) 
    n_aa::Int64 = length(individual_genotypes) - n_AA - n_Aa
    
    mean_offspring_per_genotype::Vector{Float64} = 2*exp(reproductive_rate*(1 - number_of_inds_in_patch/carring_capacity)) .* fittnesses_by_genotype

    relized_number_of_offspring_by_genotype::Vector{Int64} = [rand(Poisson((mean_offspring_per_genotype[1]  * n_AA))),
                                                              rand(Poisson((mean_offspring_per_genotype[2]  * n_Aa))),
                                                              rand(Poisson((mean_offspring_per_genotype[3]  * n_aa)))]
    
    relized_number_of_offspring::Int64 = sum(relized_number_of_offspring_by_genotype)
    number_of_allele_A::Int64 = relized_number_of_offspring_by_genotype[1] + sum(bitrand(relized_number_of_offspring_by_genotype[2]))

    return shuffle!(vcat(repeat([0],number_of_allele_A), repeat([1],relized_number_of_offspring-number_of_allele_A)))
end 

function poisson_gamete_production2(individuals::Vector{Tuple{Bool,Bool}}, number_of_inds_in_patch::Int64, location::Int64, fittnesses_by_genotype_first_half::Vector{Float64}, number_of_patches::Int64, carring_capacity::Int64,reproductive_rate::Float64)::BitVector
    #this function produces the gametes of the next generation
    #the "middle" of the habiat is at number_of_patches ÷ 2 + 0.5
    fittnesses_by_genotype::Vector{Float64} = (location <= number_of_patches ÷ 2) ? fittnesses_by_genotype_first_half :
                                                                                    reverse(fittnesses_by_genotype_first_half)
    individual_genotypes::Vector{Int64} = sum.(individuals) 
    n_AA::Int64 = sum(individual_genotypes .== 0) 
    n_Aa::Int64 = sum(individual_genotypes .== 1) 
    n_aa::Int64 = length(individual_genotypes) - n_AA - n_Aa
    
    mean_offspring_per_genotype::Vector{Float64} = 2*exp(reproductive_rate*(1 - number_of_inds_in_patch/carring_capacity)) .* fittnesses_by_genotype

    relized_number_of_offspring_by_genotype::Vector{Int64} = [rand(Poisson((mean_offspring_per_genotype[1]  * n_AA))),
                                                              rand(Poisson((mean_offspring_per_genotype[2]  * n_Aa))),
                                                              rand(Poisson((mean_offspring_per_genotype[3]  * n_aa)))]
    
    relized_number_of_offspring::Int64 = sum(relized_number_of_offspring_by_genotype)
    number_of_allele_A::Int64 = relized_number_of_offspring_by_genotype[1] + sum(bitrand(relized_number_of_offspring_by_genotype[2]))
    #number_of_allele_A::Int64 = relized_number_of_offspring_by_genotype[1] + relized_number_of_offspring_by_genotype[2]/2
    #p::Float64 =(relized_number_of_offspring_by_genotype[1] + relized_number_of_offspring_by_genotype[2]/2)/ relized_number_of_offspring
    #return rand(relized_number_of_offspring) .>p
    

    return shuffle!(vcat(repeat([0],number_of_allele_A), repeat([1],relized_number_of_offspring-number_of_allele_A)))
end 
pop = create_population(100,100)
patch = pop[100]
inds = patch.males
selection_coefficient = 0.1
fittnesses_by_genotype_first_half = [1, 1-selection_coefficient, 1-2*selection_coefficient]

histogram([mean(poisson_gamete_production(inds, 100, 1, fittnesses_by_genotype_first_half, 100,100,1.2)) for i=1:1000])


A = [mean(poisson_gamete_production1(inds, 100, 1, fittnesses_by_genotype_first_half, 100,100,1.2)) for i=1:1000]
A1 = [length(poisson_gamete_production1(inds, 100, 1, fittnesses_by_genotype_first_half, 100,100,1.2)) for i=1:1000]

mean(A1)
std(A)
A = [mean(poisson_gamete_production1(inds, 100, 1, fittnesses_by_genotype_first_half, 100,100,1.2)) for i=1:1000]
B = [mean(poisson_gamete_production2(inds, 100, 1, fittnesses_by_genotype_first_half, 100,100,1.2)) for i=1:1000]

B1 = [length(poisson_gamete_production2(inds, 100, 1, fittnesses_by_genotype_first_half, 100,100,1.2)) for i=1:1000]
mean(B1)
std(B)

histogram([A,B], bins = 20, opacity = 0.5)
histogram!(A)







function poisson_mating(patch::Patch, fittnesses_by_genotype_first_half::Vector{Float64}, number_of_patches::Int64, carring_capacity::Int64,reproductive_rate::Float64)::Vector{Tuple{Bool,Bool}}
    number_of_inds_in_patch::Int64 = length(patch.males) + length(patch.females)
    male_gametes::BitVector  = poisson_gamete_production(patch.males, number_of_inds_in_patch, patch.location, fittnesses_by_genotype_first_half, number_of_patches, carring_capacity, reproductive_rate)
    female_gametes::BitVector = poisson_gamete_production(patch.females, number_of_inds_in_patch, patch.location, fittnesses_by_genotype_first_half, number_of_patches, carring_capacity, reproductive_rate)
    number_offspring::Int64 = min(length(male_gametes),length(female_gametes))
    return tuple.(male_gametes[1:number_offspring],female_gametes[1:number_offspring])
end

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


using BenchmarkTools
individuals = patch.males
@benchmark poisson_gamete_production(inds, 200, 1, fittnesses_by_genotype_first_half,100,200,1.2)
pop = create_population(100,100)
patch = pop[100]
selection_coefficient = 0.1
fittnesses_by_genotype_first_half = [1, 1-selection_coefficient, 1-2*selection_coefficient]
mating_function = ((patch)->poisson_mating(patch, fittnesses_by_genotype_first_half, 100, 200, 1.2))

sizes = []
freqs = []
for i=1:1000
    push!(sizes, length(patch.males) + length(patch.females))
    push!(freqs, (mean(mean.(patch.males)) + mean(mean.(patch.females)))/2 )
    patch = local_mating(patch, mating_function)
end
patch.males

using Plots 

plot(freqs)

pop = create_population(100,100)
patch = pop[1]
selection_coefficient = 0.1
fittnesses_by_genotype_first_half = [1, 1-selection_coefficient, 1-2*selection_coefficient]
number_of_patches = 100
location1 = 50
individuals = patch.males
number_of_inds_in_patch = length(patch.males) + length(patch.females)
#this function produces the gametes of the next generation
#the "middle" of the habiat is at number_of_patches ÷ 2 + 0.5
fittnesses_by_genotype::Vector{Float64} = (location1 <= number_of_patches ÷ 2) ? fittnesses_by_genotype_first_half :
                                                                                reverse(fittnesses_by_genotype_first_half)
individual_genotypes::Vector{Int64} = sum.(individuals) 
n_AA::Int64 = sum(individual_genotypes .== 0) 
n_Aa::Int64 = sum(individual_genotypes .== 1) 
n_aa::Int64 = length(individual_genotypes) - n_AA - n_Aa

mean_offspring_per_genotype::Vector{Float64} = 2*exp(reproductive_rate*(1 - number_of_inds_in_patch/carring_capacity)) .* fittnesses_by_genotype

relized_number_of_offspring_by_genotype::Vector{Int64} = [rand(Poisson((mean_offspring_per_genotype[1]  * n_AA))),
                                                            rand(Poisson((mean_offspring_per_genotype[2]  * n_Aa))),
                                                            rand(Poisson((mean_offspring_per_genotype[3]  * n_aa)))]

relized_number_of_offspring::Int64 = sum(relized_number_of_offspring_by_genotype)
number_of_allele_A::Int64 = relized_number_of_offspring_by_genotype[1] + sum(bitrand(relized_number_of_offspring_by_genotype[2]))



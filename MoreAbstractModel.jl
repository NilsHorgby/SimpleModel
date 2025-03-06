using Distributions
using StatsBase

dispersal_distance::Float64 = 1.0
number_of_patches::Int64 = 100

function create_subpopulation(individuals_per_patch::Int64, p::Float64 = 0.5)::Vector{Int64}
    N_A::Int64 = rand(Binomial(2*individuals_per_patch,p))
    p_A::Float64 = N_A/(2*individuals_per_patch)
    N_a::Int64 = 2*individuals_per_patch - N_A
    p_a::Float64 = N_a/(2*individuals_per_patch)
    N_AA::Int64 = rand(Binomial(N_A,p_A^2))
    N_Aa::Int64 = rand(Binomial(N_A-2*N_AA,2*p_A*p_a))
    N_aa::Int64 = individuals_per_patch - N_AA - N_Aa
    return [N_AA, N_Aa, N_aa]
end


function create_population(number_of_patches::Int64, individuals_per_patch::Int64, p::Float64 = 0.5)::Array{Int64}
    return reduce(hcat,[create_subpopulation(individuals_per_patch,p) for i in 1:number_of_patches])
end

males::Array{Int64} = create_population(100,100)
females::Array{Int64} = create_population(100,100)
males

function discrete_distrubution_PDF(dispersal_distance::Real, distance::Int64)::Float64
    return .5 * (erf((distance + 1/2)/(sqrt(2)*dispersal_distance)) - erf((distance - 1/2)/(sqrt(2)*dispersal_distance)))
end

function discrete_distrubution_RNG(number_of_inds::Int64, dispersal_distance::Real, location::Int64)::Vector{Int64}
    try
        return round.(Int,rand(Normal(location,dispersal_distance), number_of_inds))
    catch
        println("Error, $number_of_inds, $dispersal_distance, $location")
        return Vector{Int64}(undef, number_of_inds)
    end
end


function migrate_from_to(num_one_genotype::Int64, location::Int64)::Vector{Int64}
    next_patch = discrete_distrubution_RNG(num_one_genotype, dispersal_distance, location)
    next_patch[next_patch .<= 0] .= 1
    next_patch[next_patch .> number_of_patches] .= number_of_patches
    return counts(next_patch, 1:number_of_patches)
end

function migrate_from_to(args::Tuple{Int64,Int64})::Vector{Int64}
    num_one_genotype::Int64, location::Int64 = args[1], args[2]
    next_patch = discrete_distrubution_RNG(num_one_genotype, dispersal_distance, location)
    next_patch[next_patch .<= 0] .= 1
    next_patch[next_patch .> number_of_patches] .= number_of_patches
    return counts(next_patch, 1:number_of_patches)
end
#this one fucked
function migrate_one_genotype(individuals::Vector{Int64})::Vector{Int64}
    return mapreduce(migrate_from_to,((x,y)->x[1].+y),zip(individuals,1:number_of_patches))
end
using Plots
inds = zeros(Int64,100)
inds[50] = 100
sum(inds)
p = plot()
for i in 1:100
    inds = migrate_one_genotype(inds)
    plot!(p,inds)
end
xlims!(40,60)
display(p)
sum(inds)
sum.(migrate_from_to.(inds,1:100))
sum.(migrate_from_to.(inds,1:100))





function migrate_all_genotypes(individuals::Array{Int64})
    return mapreduce(n -> migrate_one_genotype(individuals[n,:]),hcat,1:3)
end

function migration(males::Array{Int64}, females::Array{Int64})
    return sum.(migrate_from_to.(males,50)), sum.(migrate_from_to.(females,50))
end
quit = exit
using Distributions
using Distributed
using Plots
addprocs(7)
cd("/home/gadus/programing-projects/julia-projects/SimpleModel")

@everywhere include("Initialisation.jl")
@everywhere include("Mating.jl")
@everywhere include("Migration.jl")
@everywhere using .Initialisation, .Mating, .Migration

@everywhere using Statistics: mean, std


@everywhere function run_simulation(selection_cofficient::Float64, dispersal_distance::Float64, number_of_patches::Int64, individuals_per_patch::Int64)
    population = create_population(number_of_patches, individuals_per_patch)
    for n in 1:100
        population = global_mating(population, selection_cofficient, number_of_patches, individuals_per_patch) |> (x -> migration(x,dispersal_distance,number_of_patches))
    end
    return population
end


@everywhere function calculate_allele_frequency(patch:: Patch)
    return mean((Individual -> mean(Individual.genome)).(patch.females))/2 + mean((Individual -> mean(Individual.genome)).(patch.males))/2
end

@everywhere function get_allele_frequencies(selection_cofficient::Float64, dispersal_distance::Float64,  number_of_patches::Int64, individuals_per_patch::Int64)::Vector{Float64}
    population = run_simulation(selection_cofficient, dispersal_distance, number_of_patches, individuals_per_patch)
    #This line is awful, git gud ig
    calculate_allele_frequency = patch -> mean((Individual -> mean(Individual.genome)).(patch.females))/2 + mean((Individual -> mean(Individual.genome)).(patch.males))/2
    return calculate_allele_frequency.(population)
end


function predicted_allele_frequency(patch_number, selection_cofficient::Float64, dispersal_distance::Float64, number_of_patches::Int64)::Float64
    x = (patch_number+.5) - number_of_patches/2
    if x >= 0
        return -1/2 + 3/2 * tanh(sqrt(selection_cofficient/(2*dispersal_distance^2))*x+ atanh(sqrt(2/3)))^2
    else
        return 3/2 - 3/2 * tanh(-sqrt(selection_cofficient/(2*dispersal_distance^2))*x + atanh(sqrt(2/3)))^2
    end
end

@everywhere function get_sum_pq(selection_cofficient::Float64, dispersal_distance::Float64, n_runs::Int64)::Tuple{Tuple{Float64,Float64},Float64}
    sum_pq = []
    mean_p = zeros(number_of_patches)
    for n in 1:n_runs
        population = run_simulation(selection_cofficient, dispersal_distance, 100, 100)
        p = calculate_allele_frequency.(population)
        q = 1 .- p
        mean_p .+= p / n_runs
        push!(sum_pq,sum(p .* q))
    end
    return (4*mean(sum_pq),std(4*sum_pq)), 4*sum(mean_p .* (1 .- mean_p))
end

function calculate_cline_width(selection_cofficient::Float64, dispersal_distance::Float64)
    return dispersal_distance/sqrt(selection_cofficient/3)
end

function calculate_cline_width(input::Tuple{Float64,Float64})
    selection_cofficient, dispersal_distance = input
    return dispersal_distance/sqrt(selection_cofficient/3)
end


function get_allele_frequencies(population::Vector{Patch})::Vector{Float64}
    #This line is awful, git gud ig
    calculate_allele_frequency = patch -> mean((Individual -> mean(Individual.genome)).(patch.females))/2 + mean((Individual -> mean(Individual.genome)).(patch.males))/2
    return calculate_allele_frequency.(population)
end


function plot_populations(populations::Vector{Vector{Patch}},s::Float64,l::Float64,number_of_patches::Int64,individuals_per_patch::Int64, n_runs::Int64)
    frequencies = get_allele_frequencies.(populations)
    print(number_of_patches)
    p = plot()
    plot!(p, frequencies,linestyle = :dashdot, color = "gray", label = "")
    plot!(p, mean.(collect(eachrow(reduce(hcat, frequencies)))), 
        title = "Allele frequency in Simple Model, mean of $n_runs runs\n s = $s, l = $l, N = $individuals_per_patch",
        xlabel = "Patch",
        ylabel = "Allele frequency",
        line = (3,:black,:dashdot),
        label = "Mean Observed Allele Frequencies")
    actuall_l = std(round.(rand(Normal(0,l),100000)))
    predicted_frequencies = [predicted_allele_frequency(i,s,actuall_l,number_of_patches) for i in range(1,number_of_patches,length = 1000)]
    plot!(p, range(1,Int(number_of_patches),length = 1000),predicted_frequencies,
        label = "Predicted allele frequencies",
        line=(3,:red))
    return p
end

@everywhere function run_simulation_sum_pqs(selection_cofficient::Float64, dispersal_distance::Float64, number_of_patches::Int64, individuals_per_patch::Int64) :: Tuple{Vector{Float64},Vector{Patch}}
    population = create_population(number_of_patches, individuals_per_patch)
    sum_pq = []
    for n in 1:10000
        population = global_mating(population,
                                   selection_cofficient,
                                   number_of_patches,
                                   individuals_per_patch) |> (x -> migration(x,dispersal_distance,number_of_patches))
        p = calculate_allele_frequency.(population)
        if n%10==0
            push!(sum_pq, 4.0 * sum(p .* (1 .- p)))
        end
    end
    return sum_pq, population
end
p = [1,2,3]

x = pmap(x->run_simulation_sum_pqs(0.05,0.3,100,1000),1:30)

sum_pqs = (x -> x[1]).(x)
population = (x -> x[2]).(x)

p = plot_populations(population,0.05,0.3,100,1000,30)

plot!([0,100],[0.2,0.2], label="")
plot!([0,100],[0.8,0.8], label="")
xlims!(p, 45, 55)


traveling_mean = map(x -> [mean(x[900:i]) for i in 900:1000], sum_pqs)


p = plot()
generations = 9000:10:10000


plot!(p,generations,3/4*(mean(traveling_mean)+std(traveling_mean)),
      label = "Mean + std",line =  (:gray))
plot!(p,generations,3/4*mean(traveling_mean), label = "Mean Cline Width")
plot!(p,generations,3/4*(mean(traveling_mean)-std(traveling_mean)), label = "Mean - std",line =  (:gray))

plot!(p, [9000,10000], [calculate_cline_width(0.05,0.3),calculate_cline_width(0.05,0.3)], label = "Predicted Cline Width", line= (3,:black,:dashdot))

savefig(p, "figures/cline-width-traveling-mean-log-scale-larger-pop.png")


calculate_cline_width(0.05,0.3)
3/4*mean(traveling_mean)[end]


using StatsBase: sample, Weights

mean(sample([0,1], Weights([2,3]), 100000))
3/5
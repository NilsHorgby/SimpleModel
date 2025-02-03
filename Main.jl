quit = exit
using Distributed
addprocs(4)

@everywhere include("Initialisation.jl")
@everywhere include("Mating.jl")
@everywhere include("Migration.jl")
@everywhere using .Initialisation, .Mating, .Migration
@everywhere using Statistics: mean, std

@everywhere function run_simulation(selection_cofficient::Float64)
    population = Initialisation.create_population()
    for n in 1:100
        population = Mating.global_mating(population, selection_cofficient) |> Migration.migration
    end
    return population
end
population = create_population()

patch = population[1]

@everywhere function calculate_allele_frequency(patch:: Patch)
    return mean((Individual -> mean(Individual.genome)).(patch.females))/2 + mean((Individual -> mean(Individual.genome)).(patch.males))/2
end

@everywhere function get_allele_frequencies(selection_cofficient::Float64)::Vector{Float64}
    population = run_simulation(selection_cofficient)
    #This line is awful, git gud ig
    calculate_allele_frequency = patch -> mean((Individual -> mean(Individual.genome)).(patch.females))/2 + mean((Individual -> mean(Individual.genome)).(patch.males))/2
    return calculate_allele_frequency.(population)
end


function predicted_allele_frequency(patch_number, selection_cofficient::Float64)::Float64

    x = patch_number - number_of_patches/2
    if x >= 0
        return -1/2 + 3/2 * tanh(sqrt(selection_cofficient)*x/(2*1.5) + atanh(sqrt(2/3)))^2
    else
        return 3/2 - 3/2 * tanh(-sqrt(selection_cofficient)*x/(2*1.5) + atanh(sqrt(2/3)))^2
    end
end

frequencies = pmap(x -> get_allele_frequencies(0.02), 1:100)

using Plots


plot(frequencies,linestyle = :dashdot, color = "gray", label = "")
plot!(mean.(collect(eachrow(reduce(hcat, frequencies)))), 
    title = "Allele frequency in Simple Model,\n mean of 100 runs, s = 0.02, l = 1.5",
    xlabel = "Patch",
    ylabel = "Allele frequency",
    line = (3,:black,:dashdot),
    label = "Mean Observed Allele Frequencies")

predicted_frequencies = [predicted_allele_frequency(i,0.2) for i in range(1,Int(number_of_patches),length = 1000)]
plot!(range(1,Int(number_of_patches),length = 1000),predicted_frequencies,
    label = "Predicted allele frequencies",
    line=(3,:red))
savefig("~/programing-projects/julia-projects/SimpleModel/BasicModel_with_pred_s0-2_l1-5.png")
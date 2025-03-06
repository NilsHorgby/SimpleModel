using Distributions
cd("/home/gadus/programing-projects/julia-projects/SimpleModel")
include("Initialisation.jl")
include("Mating.jl")
include("Migration.jl")
include("Utils.jl")
using .Initialisation, .Mating, .Migration, .Utils
using Statistics: median, mean, std
using Base.Threads
using JLD2






function calculate_allele_frequency(patch::Patch)::Vector{Int64}
    number_of_AA::Int64 = sum(sum.(patch.males )  .==2) +
                          sum(sum.(patch.females ).==2)
    number_of_Aa::Int64 = sum(sum.(patch.males)   .==1) +
                          sum(sum.(patch.females) .==1)
    number_of_aa::Int64 = sum(sum.(patch.males)   .==0) +
                          sum(sum.(patch.females) .==0)
    return [number_of_AA,number_of_Aa,number_of_aa]
end

function run_simulation(dispersal_distance::Real, gamete_producing_function::Function, generations::Int64, number_of_patches::Int64, individuals_per_patch::Int64, sampeling_frequency::Int64)
    population = create_population(number_of_patches, individuals_per_patch,0.5)
    allele_frequencies_over_time:: Array{Float64} = zeros(number_of_patches,3,ceil(Int,2*generations/sampeling_frequency))
    
    for gen in 1:generations
        population = global_mating(population, gamete_producing_function)
        if gen % sampeling_frequency == 0
            #allele_frequencies_over_time[:,:,floor(Int,2*gen/sampeling_frequency)-1] = reduce(vcat,calculate_allele_frequency.(population)')
            population = migration(population, dispersal_distance)
            allele_frequencies_over_time[:,:,floor(Int,2*gen/sampeling_frequency)] = reduce(vcat,calculate_allele_frequency.(population)')
        else
            population = migration(population, dispersal_distance)
        end
    end
    return allele_frequencies_over_time
end

runs::Int64 = 10
generations::Int64 = 10000
number_of_patches::Int64 = 100
individuals_per_patch::Int64 = 100
carring_capacity::Int64 = 200
selection_coefficients::Vector{Float64} = [0.1,0.2,0.3]
cline_widths::Vector{Int64} = [10,15,20]
reproductive_rates::Vector{Float64} = [1.1,1.25]
sampeling_frequency::Int64 = 2

dispersal_distance::Float64 = sqrt(cline_widths[1] ^ 2 * selection_coefficients[1]/3)
fittnesses_by_genotype_first_half = (1.0,1-selection_coefficients[1],1-2*selection_coefficients[1])
curried_mating = curry_gamete_production(poisson_mating, fittnesses_by_genotype_first_half, number_of_patches, carring_capacity, 1.2)
curried_mating = curry_gamete_production(multinomial_mating, fittnesses_by_genotype_first_half, number_of_patches, individuals_per_patch)

run_simulation(dispersal_distance,curried_mating,generations,number_of_patches,individuals_per_patch, sampeling_frequency)


full_model_consts = quote
    const runs::Int64 = 30
    const generations::Int64 = 10_000
    const number_of_patches::Int64 = 300
    const individuals_per_patch::Int64 = 100
    const carring_capacity::Int64 = 100
    const selection_coefficients::Vector{Float64} = [0.1,0.2,0.3]
    const cline_widths::Vector{Int64} = [10,15,20]
    const reproductive_rates::Vector{Float64} = [1.1,1.25]
    const sampeling_frequency::Int64 = 20
end

test_model_consts = quote
    const runs::Int64 = 10
    const generations::Int64 = 10
    const number_of_patches::Int64 = 10
    const individuals_per_patch::Int64 = 100
    const carring_capacity::Int64 = 200
    const selection_coefficients::Vector{Float64} = [0.1,0.2,0.3]
    const cline_widths::Vector{Int64} = [10,15,20]
    const reproductive_rates::Vector{Float64} = [1.1,1.25]
    const sampeling_frequency::Int64 = 2
end
if ARGS[1] == "full"
    eval(full_model_consts)
elseif ARGS[1] == "test"
    eval(test_model_consts)
else
    print("""Did not recoginice argument, please selected argument in {"full", "test"} """)
end
global n_left::Int64 = length(selection_coefficients)*length(cline_widths)*length(reproductive_rates)

function main(dir::String,mating_type::String) 
    t00::Float64 = time()
    times::Vector{Float64} = []
    for cline_width in cline_widths; for selection_coefficient in selection_coefficients; for reproductive_rate in reproductive_rates
        if ~isfile("$dir/allele_freqs_by_time_s_$(selection_coefficient)_cw_$(cline_width)_r_$reproductive_rate.jld2")
            println("Doing s=$(selection_coefficient) cw=$(cline_width) r=$reproductive_rate")
            t0 = time()
            fittnesses_by_genotype_first_half::Tuple{Float64,Float64,Float64} = (1,1-selection_coefficient,1-2*selection_coefficient)
            dispersal_distance::Float64 = sqrt(cline_width ^ 2 * selection_coefficient/3)
            

            curried_mating::Function = (mating_type == "Poisson") ? curry_gamete_production(poisson_mating, fittnesses_by_genotype_first_half, number_of_patches, carring_capacity, reproductive_rate) :
                                                                    curry_gamete_production(multinomial_mating, fittnesses_by_genotype_first_half, number_of_patches, individuals_per_patch)
 
            allele_frequencies_over_time::Array{Float64} = zeros(runs,number_of_patches,3,ceil(Int,2*generations/sampeling_frequency))

            @threads for run in 1:runs
                allele_frequencies_over_time[run,:,:,:] = run_simulation(dispersal_distance,curried_mating,generations,number_of_patches,individuals_per_patch, sampeling_frequency)
            end
            if ARGS[1] == "full"
                save_object("$dir/allele_freqs_by_time_s_$(selection_coefficient)_cw_$(cline_width)_r_$reproductive_rate.jld2",allele_frequencies_over_time)
            else
                println("""Did not save, input "full" if this is desiered """)
            end
            push!(times, abs(time() - t0))
            median_time_left::String = format_time(median(times) * n_left)
            proc_time_left::String = format_time(1.96 * std(times)*n_left/sqrt(length(times)))
            println("Done with s=$(selection_coefficient) cw=$(cline_width) r=$reproductive_rate, Took $(times[end]), Time left $(median_time_left) +- $(proc_time_left)")  
        else
            println("Skipping s=$(selection_coefficient) cw=$(cline_width) r=$reproductive_rate")
        end
    end;end;end
end


main("results-poisson-model","Poisson")
main("results-binomial-model", "Binomial")


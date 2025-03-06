using Distributions
cd("/home/gadus/programing-projects/julia-projects/SimpleModel")
include("Initialisation.jl")
include("Mating.jl")
include("Migration.jl")
using .Initialisation, .Mating, .Migration
using Statistics: median, mean, std
using Base.Threads


full_model_consts = quote
    const number_of_patches::Int64 = 300
    const individuals_per_patch::Int64 = 100
    const carring_capacity::Int64 = 100
    const reproductive_rate::Float64 = 1.2 #look at effect of large rm
    const runs::Int64 = 100
    #width of cline should be 10 - 30 

    const number_of_patches::Int64 = 100
    const untracked_generations::Int64 = 9000
    const tracked_generations::Int64 = 1000
end

test_model_consts = quote
    const number_of_patches::Int64 = 100
    const individuals_per_patch::Int64 = 300
    const carring_capacity::Int64 = 300
    const reproductive_rate::Float64 = 1.2
    const runs::Int64 = 10

    const number_of_patches::Int64 = 100
    const untracked_generations::Int64 = 90
    const tracked_generations::Int64 = 10
endpoisson_model_consts
if ~isinteractive() 
if ARGS[1] == "full"
    eval(full_model_consts)
elseif ARGS[1] == "test"
    eval(test_model_consts)
else
    println("Please provide either 'full' or 'test' as the first argument, to specify the scale of the model to run")
    exit()
end

poisson_model_consts = quote
    #const gamete_producing_function = poisson_gamete_production
    const dir = "results-poisson-model"
end

binomial_model_consts = quote
    #const gamete_producing_function = binomial_gamete_production
    const dir = "results-binomial-model"
end

if ARGS[2] == "poisson"
    eval(poisson_model_consts)
elseif ARGS[2] == "binomial"
    eval(binomial_model_consts)
else
    println("Please provide either 'poisson' or 'binomial' as the second argument, to specify the gamete production model to use")
    exit()
end
end



function curry_gamete_production(f::Function, args...)::Function
    return (inds, location)->f(inds,location,args...)
end

function run_simulation(dispersal_distance::Real, gamete_producing_function::Function)
    population = create_population(number_of_patches, individuals_per_patch)
    population_size_over_time::Vector{Vector{Int64}} = Vector(undef,ceil(Int,untracked_generations/10))
    
    for gen in 1:untracked_generations
        population = global_mating(population, gamete_producing_function)
        population = migration(population, dispersal_distance)
        if gen % 10 == 0
            population_size_over_time[floor(Int,gen/10)] = [length(patch.males) + length(patch.females) for patch in population]
        end
    end
    return population,population_size_over_time
end


function calculate_allele_frequency(patch::Patch)::Float64
    return mean(mean.(patch.males))/2 + mean(mean.(patch.females))/2
end

function run_simulation_and_track_pq(population::Vector{Patch}, population_size_over_time::Vector{Vector{Int64}}, dispersal_distance::Float64, gamete_producing_function::Function)::Tuple{Array{Float64},Vector{Int64}}
    allele_frequencies_over_time:: Array{Float64} = zeros(number_of_patches,tracked_generations)
    
    for gen in 1:tracked_generations
        population = global_mating(population, gamete_producing_function)
        population = migration(population, dispersal_distance)
        allele_frequencies_over_time[:,gen] = calculate_allele_frequency.(population)
        if gen % 10 == 0
            push!(population_size_over_time, [length(patch.males) + length(patch.females) for patch in population])
        end
    end
    return allele_frequencies_over_time, population_size_over_time
end




function run_simulation_and_stability_analysis(dispersal_distance::Float64, gamete_producing_function::Function)
    populations = Vector{Vector{Patch}}(undef, runs)
    population_sizes_over_time = Vector{Vector{Int64}}(undef, runs)
    allele_frequencies_over_time::Array{Float64} = zeros(runs, number_of_patches, tracked_generations)

    @threads for i in 1:runs
        populations[i], population_sizes_over_time[i] = run_simulation(dispersal_distance, gamete_producing_function)
        allele_frequencies_over_time[i,:,:], population_sizes_over_time[i] = run_simulation_and_track_pq(populations[i],
                                                                                                         population_sizes_over_time[i],
                                                                                                         dispersal_distance,
                                                                                                         gamete_producing_function)
    end
    return allele_frequencies_over_time, population_sizes_over_time
end

using JLD2

function format_time(time_s::Float64)::String
    if isnan(time_s)
        return "NaN"
    else
        hours = floor(Int,time_s/3600)
        minutes = floor(Int,(time_s - hours*3600)/60)
        seconds = round(Int,time_s%60)
        return "$hours:$minutes:$seconds"
    end
end
const s::Vector{Float64}= [0.05,0.1,0.2,0.3,0.4,0.5]
const l::Vector{Float64}= [0.3,0.6,1,1.5,2,3]

global n_left::Int64 = length(s)*length(l)

#main function: Method of gamete production, directory to save files
function main(dir::String) 
    t00::Float64 = time()
    times::Vector{Float64} = []
    
    for si in s
        for li in l
            t0::Float64 = time()
            global n_left -= 1
            if ~isfile("$dir/frequencies_s_$si"*"_l_$li.jld2")
                println("Doing s = $si l = $li")
                if ARGS[2] == "poisson"
                    curried_gamete_producing_function = curry_gamete_production(poisson_gamete_production, si, number_of_patches, carring_capacity, reproductive_rate)
                elseif ARGS[2] == "binomial"
                    curried_gamete_producing_function = curry_gamete_production(binomial_gamete_production, si, number_of_patches, individuals_per_patch)
                else
                    println("This should not happen, ARGS check failed")
                    exit()
                end
                
                allele_frequencies_over_time, population_sizes_over_time = run_simulation_and_stability_analysis(li, curried_gamete_producing_function)
                
                if ARGS[1] == "full"
                save_object("$dir/frequencies_s_$si"*"_l_$li.jld2",allele_frequencies_over_time)
                save_object("$dir/sizes_s_$si"*"_l_$li.jld2",population_sizes_over_time)
                end
                push!(times, abs(time() - t0))
                median_time_left::String = format_time(median(times) * n_left)
                proc_time_left::String = format_time(1.96 * std(times)*n_left/sqrt(length(times)))
                println("Done with s = $si l = $li, took $(format_time(times[end])), time left: $median_time_left ± $proc_time_left")
            else
                println("Skipping populations_s_$si"*"_l_$li.jld2")
            end
        end
    end
    total_time::Float64 = time()-t00
    if ARGS[1] == "full"
    save_object("$dir/times_taken.jld2",times)
    save_object("$dir/parameters.jld2",Dict("s"=>s,"l"=>l))
    save_object("$dir/total_time.jld2",total_time)
    end
    println("Total time: $(format_time(total_time))")
end
main(dir)
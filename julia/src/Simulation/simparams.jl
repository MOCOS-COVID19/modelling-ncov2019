using Random
using Distributions
rng=MersenneTwister(0)
using CSV


struct RunParams
  seed::Int

  constant_kernel_param::Float64
  household_kernel_param::Float64

  backward_tracking_prob::Float32
  backward_detection_delay::TimeDiff

  forward_tracking_prob::Float32
  forward_detection_delay::TimeDiff

  quarantine_length::Float32
  testing_time::TimeDiff
end

struct PopulationParams
  household_ptrs::Vector{Tuple{UInt32,UInt32}}  # (i1,i2) where i1 and i2 are the indices of first and last member of the household
  
end

struct Progression
    severity::Severity
    # times given with respect to the infection time
    incubation_time::Union{Missing,TimeDiff}
    mild_symptoms_time::Union{Missing,TimeDiff}
    severe_symptoms_time::Union{Missing,TimeDiff}
    #critical_symptoms_time::Float32
    recovery_time::Union{Missing,TimeDiff}
    death_time::Union{Missing, TimeDiff}
    #Progression(severity::Severity, incubation_time::Real, mild_time::Real, severe_time::Real, recovery_time) = incubation < mild_time < severe_time < recovery_time
end


function create_dist_severity(age_induced_fatality_rates) 
  dist =[] #Array{Float32}
  for mod_p ∈ age_induced_fatality_rates[:, 3]
      age_dep = mod_p/death_probability.Critical.p
      hist = def_dist_severity.p ./ (1 - def_dist_severity.p[4]) .* (1 - age_dep )
      hist[4] = age_dep
      push!(dist, Categorical(hist))
  end
  return dist
end

#= # test 
dist_severity = create_dist_severity(age_induced_fatality_rates)
for i ∈  dist_severity
    println(i.p)
end
output:
[0.006270210077735186, 0.8454333254812944, 0.1442148317879093, 0.004081632653061225]
[0.006270210077735186, 0.8454333254812944, 0.1442148317879093, 0.004081632653061225]
[0.006244512495449387, 0.8419684348030924, 0.1436237873953359, 0.00816326530612245]
[0.006128873375163287, 0.8263764267511833, 0.14096408762875562, 0.026530612244897958]
[0.005833351178876588, 0.78653018395186, 0.13416707711416154, 0.07346938775510203]
[0.005268004368588989, 0.7103025890314154, 0.12116410047754676, 0.163265306122449]
[0.00439428657087179, 0.5924963059725464, 0.10106859113005119, 0.3020408163265306
# VS oryginal from python
#[0.00627021 0.84543333 0.14421483 0.00408163]
#[0.00627021 0.84543333 0.14421483 0.00408163]
#[0.00624451 0.84196843 0.14362379 0.00816327]
#[0.00612887 0.82637643 0.14096409 0.02653061]
#[0.00583335 0.78653018 0.13416708 0.07346939]
#[0.005268   0.71030259 0.1211641  0.16326531]
#[0.00439429 0.59249631 0.10106859 0.30204082]
 =#


function get_severity(dist_severity)
  ☠ = false
  severity = rand(rng, dist_severity) |> Severity
  if severity == Critical
      dist_age_fatality_rate = Bernoulli(death_probability.Critical.p)
      ☠ = rand(rng, dist_age_fatality_rate)
  end
  return severity, ☠
end


function sample_progression(rng::AbstractRNG, dist_severity, dist_symptom_onset,
  dist_stay_home, dist_hospitalization, dist_☠, age)

  t₀   = rand(rng, dist_symptom_onset) # + t₋₁
  t₁   = t₀ + rand(rng, dist_stay_home)
  t₂   = t₀ + rand(rng, dist_hospitalization)
  t☠ = NaN

  #severity = rand(rng, dist_severity) |> Severity
  severity, ☠ = get_severity(dist_severity)
  mild_symptoms_time = Inf
  severe_symptoms_time = Inf
  recovery_time = Inf

  incubation_time = t₀

  if (severity==Mild) || (severity==Asymptomatic)
    mild_symptoms_time = t₁
    severe_symptoms_time = NaN
    recovery_time = t₀ + 2*7
  end

  if (severity==Severe) || (severity==Critical)
    if (t₁ <= t₂)
      mild_symptoms_time = t₁
    else
      mild_symptoms_time = NaN
    end
    severe_symptoms_time = t₂
    recovery_time = t₀ + 6*7
  end

  if ☠ == true
    t☠  = t₀ + rand(rng, dist_☠)        
  end

  #println(t₀,"\t", t₁,"\t", t₂,"\t", t☠, "\t", recovery_time )

  Progression(
    severity,
    incubation_time,
    mild_symptoms_time,
    severe_symptoms_time,
    recovery_time,
    t☠
  )
end


#### How to use it! ##############


#= 
range=1000000
age = rand(81:81,range)

#make it once!
dist_severity = create_dist_severity(age_induced_fatality_rates)

a=[]
for i ∈ age 
      age_range = findfirst( age_induced_fatality_rates[:,2] .>= i )
      dist_severity_local = dist_severity[age_range]
      choosen = sample_progression(rng, dist_severity_local, dist_symptom_onset, dist_stay_home, dist_hospitalization, dist_☠)
      push!(a, [choosen.severity,choosen.death_time])
end

bool = (x->Base.isnan(x[2])).(a)
∑☠ = range - count(bool)

bool = a[1,:] .== [(Critical, )]
∑⚠ = count((x->x[1]).(a) .== Critical)

println("∑⚠ = ", ∑⚠)
println("∑☠ = ", ∑☠)

output:
∑⚠ = 301756
∑☠ = 148140 

=#

struct SimParams 
  household_ptrs::Vector{Tuple{UInt32,UInt32}}  # (i1,i2) where i1 and i2 are the indices of first and last member of the household
    
  progressions::Vector{Progression}
    
  constant_kernel_param::Float64
  household_kernel_param::Float64

  backward_tracking_prob::Float32
  backward_detection_delay::TimeDiff
  
  forward_tracking_prob::Float32
  forward_detection_delay::TimeDiff
  
  quarantine_length::Float32
  testing_time::TimeDiff
end

progressionof(params::SimParams, person_id::Integer) = params.progressions[person_id]
severityof(params::SimParams, person_id::Integer) = progressionof(params, person_id).severity
householdof(params::SimParams, person_id::Integer) = UnitRange(params.household_ptrs[person_id]...)


#### This have to be read from JSON probably#####
#### data come from python's file defaults.py  ###
death_probability=(Asymptomatic=Bernoulli(0),
                   Mild=Bernoulli(0),
                   Severe=Bernoulli(0),
                   Critical=Bernoulli(.49))
def_dist_severity = Categorical([6/1000, 809/1000, 138/1000, 47/1000])
age_induced_fatality_rates = [ 0 20  0.002;
                              20 40  0.002; 
                              40 50  0.004;
                              50 60  0.013;
                              60 70  0.036;
                              70 80  0.080;
                              80 200 0.148 ];



function load_params(rng=MersenneTwister(0);
        population::Union{AbstractString,DataFrame},
        kwargs...
        )
  individuals_df::DataFrame = isa(population, AbstractString) ? load_individuals(population) : population
  
  num_individuals = individuals_df |> nrow

  dist_severity = create_dist_severity(age_induced_fatality_rates)


  # files generated from the Python version of the project  
  #=
  import numpy as np
  t0LogNorm = np.load("test/models/assets/incubation_period_distribution.npy")
  t1Gamma    = np.load("test/models/assets/t1_distribution.npy")
  t2Gamma    = np.load("test/models/assets/t1_t2_distribution.npy")
  tDeathLogNorm = np.load("test/models/assets/onset_death_distribution.npy")

  np.savetxt("julia/data/t0LogNorm_incubation_period_distribution.csv", t0LogNorm ,delimiter=',')
  np.savetxt("julia/data/t1Gamma_distribution.csv", t1Gamma ,delimiter=',')
  np.savetxt("julia/data/t2Gamma_distribution.csv", t2Gamma ,delimiter=',')
  np.savetxt("julia/data/tDeathLogNorm_distribution.csv", tDeathLogNorm ,delimiter=',')
  =#


  t0LogNorm = CSV.read("julia/data/t0LogNorm_incubation_period_distribution.csv")
  t1Gamma = CSV.read("julia/data/t1Gamma_distribution.csv")
  t2Gamma = CSV.read("julia/data/t2Gamma_distribution.csv")
  t☠LogNorm = CSV.read("julia/data/tDeathLogNorm_distribution.csv")

  dist_symptom_onset = fit(LogNormal, t0LogNorm[!,1])
  dist_stay_home = fit(Gamma, t1Gamma[!,1])
  dist_hospitalization = fit(Gamma, t2Gamma[!,1])
  dist_☠ = fit(LogNormal, t☠LogNorm[!,1])

#TODO#
  progressions = individuals_df.age .|> age -> sample_progression(rng::AbstractRNG, dist_severity, dist_symptom_onset,
  dist_stay_home, dist_hospitalization, dist_☠, age)
  
  make_params(rng, individuals_df=individuals_df, progressions=progressions)
end

function make_params(rng::AbstractRNG=MersenneTwister(0);
        individuals_df::DataFrame,
        progressions::AbstractArray{Progression},

        constant_kernel_param::Float64=1.0,
        household_kernel_param::Float64=1.0,
        
        backward_tracking_prob::Float64=1.0,
        backward_detection_delay::Float64=1.0,
        
        forward_tracking_prob::Float64=1.0,
        forward_detection_delay::Float64=1.0,
        
        quarantine_length::Float64=14.0,
        testing_time::Float64=1.0
        )
  sort!(individuals_df, :household_index)

  num_individuals = individuals_df |> nrow
    
  @assert num_individuals == length(progressions)

  household_ptrs = collect( zip(groupptrs(individuals_df.household_index)...))
  
  params = SimParams(
    household_ptrs,
    progressions,        
    constant_kernel_param,   
    household_kernel_param,
    
    backward_tracking_prob,
    backward_detection_delay,
    
    forward_tracking_prob,
    forward_detection_delay,
    
    quarantine_length, # quarantine length
    testing_time # testing time
  )
  params
end
module UtilsBonito
using PiecewiseDeterministicMarkovProcesses
using ProgressMeter
using SynapseElife
const PDMP = PiecewiseDeterministicMarkovProcesses

# utils_synapse.jl — identique à l'original.
# Aucune dépendance à GLMakie : les observables Makie sont
# universelles (GLMakie et WGLMakie partagent la même API).

function _evolve_synapse(
    xc0::Vector{T},
    xd0,
    p_synapse::SynapseParams,
    events_sorted_times,
    is_pre_or_post_event,
    bap_by_epsp,
    is_glu_released,
    algos;
    progress,
    verbose = false,
    plasticity_obs,
    abstol = 1e-8,
    reltol = 1e-7,
    save_positions = (false, true),
    nu = SynapseElife.buildTransitionMatrix(),
    kwargs...,
) where {T}
    # abstol = 1e-8, reltol = 1e-7

    if save_positions isa Tuple{Bool,Bool}
        save_positionsON = save_positions
        save_positionsOFF = save_positions
        # save_positionsLAST = save_positions
    else
        save_positionsON = save_positions[1]
        save_positionsOFF = save_positions[2]
        # save_positionsLAST = save_positions
    end

    @assert eltype(is_pre_or_post_event) == Bool "Provide booleans for glutamate releases."
    @assert eltype(is_glu_released) == Bool "Provide booleans for glutamate indices."

    verbose && printstyled(color = :red, "\n+++++++++++++++++++++++++++++\n")
    verbose && printstyled(color = :red, "Synapse simulation")

    XC = SynapseElife.VectorOfArray([xc0]) # vector to hold continuous variables
    XD = SynapseElife.VectorOfArray([xd0]) # vector to hold discrete variables
    tt = [0.0]  # vector of times

    # results from the simulation of an external event
    res = PDMP.PDMPResult([0.0, 0.0], copy(XC), copy(XD))

    # we collect which external events correspond to BaPs
    events_bap = events_sorted_times[is_pre_or_post_event .== false]

    # function to simulate the synapse when Glutamate is ON
    SimGluON =
        (xc, xd, t1, t2, glu) -> SynapseElife.pdmpsynapse(
            xc,
            xd,
            t1,
            t2,
            events_bap,
            bap_by_epsp,
            glu,
            p_synapse,
            nu;
            algo = algos[1],
            save_positions = save_positionsON,
            reltol = reltol,
            abstol = abstol,
            kwargs...,
        )

    # function to simulate the synapse when Glutamate is OFF
    SimGluOFF =
        (xc, xd, t1, t2) -> SynapseElife.pdmpsynapse(
            xc,
            xd,
            t1,
            t2,
            events_bap,
            bap_by_epsp,
            zero(T),
            p_synapse,
            nu;
            algo = algos[2],
            save_positions = save_positionsOFF,
            reltol = reltol,
            abstol = abstol,
            kwargs...,
        )

    # variable to display progressbar during simulation
    # +1 for the last big till p_synapse.t_end
    pbar = Progress(length(events_sorted_times) + 1, 1)
    progress[] = 0

    # random variable for Glutamate concentration
    gluDist = SynapseElife.Gamma(1/p_synapse.glu_cv^2, p_synapse.glu_cv^2)

    # we loop over the external events, simulate them and append to res
    for (eveindex, eve) in enumerate(events_sorted_times)
        verbose && printstyled(color = :red, "\n$(eveindex) / $(length(events_sorted_times)) +++++++++++++++++++++++\n")
        if is_pre_or_post_event[eveindex] == true # it is a pre-synaptic event
            # we simulate the synapse with Glutamate OFF until event time
            # then we put  Glutamate ON for dt = p_synapse.glu_width with variable amplitude (concentration)
            verbose && SynapseElife.@printf("=> Glu Off,%4d, t ∈ [%9.4e, %9.4e]\n", eveindex, tt[end], eve)

            # simulate the event with Glutamate OFF
            res = SimGluOFF(res.xc[:, end], res.xd[:, end], tt[end], eve)
            append!(XC, res.xc)
            append!(XD, res.xd)
            append!(tt, res.time)
            gluamp = rand(gluDist)
            verbose &&
                SynapseElife.@printf("=> Glu on, %4d, t ∈ [%9.4e, %9.4e]\n", eveindex, eve, eve + p_synapse.glu_width)

            # simulate the event with Glutamate ON
            # variability here
            res = SimGluON(
                res.xc[:, end],
                res.xd[:, end],
                eve,
                eve + p_synapse.glu_width,
                ifelse(is_glu_released[eveindex], gluamp, zero(T)),
            )
            append!(XC, res.xc)
            append!(XD, res.xd)
            append!(tt, res.time)

            # update camkii observable
            #append!(plasticity_obs.val.t,res.time)
            _camk2 = SynapseElife.getCamKII(nothing, res.xc, nothing)
            #append!(plasticity_obs.val.camk2,_camk2)
            _can = SynapseElife.getCaN(nothing, res.xc, nothing)
            #append!(plasticity_obs.val.can,_can)
            _weight = res.xd[38, :] .- res.xd[37, :]
            #append!(plasticity_obs.val.weight, _weight)
            #notify(plasticity_obs)

            plasticity_obs[] = (
                t = vcat(plasticity_obs[].t, res.time),
                camk2 = vcat(plasticity_obs[].camk2, _camk2),
                can = vcat(plasticity_obs[].can, _can),
                weight = vcat(plasticity_obs[].weight, _weight),
            )

        end

        # update the progress bar
        next!(pbar; showvalues = [(:steps, eveindex), (:t, tt[end])])
        progress[] = eveindex / (length(events_sorted_times) + 1) * 100
        notify(progress)
    end

    # reaching tend: we simulate the synapse with Glutamate OFF until simulation end time required
    # by the user. In  most protocol, this is taking most of the time.
    verbose && SynapseElife.@printf("=> Reaching the end, t ∈ [%9.4e, %9.4e]\n", tt[end], p_synapse.t_end)
    res = @time SimGluOFF(res.xc[:, end], res.xd[:, end], tt[end], p_synapse.t_end)
    @info "last bit" length(res.time) tt[end] p_synapse.t_end
    append!(XC, res.xc)
    append!(XD, res.xd)
    append!(tt, res.time)

    # update camkii observable
    #append!(plasticity_obs.val.t,res.time)
    _camk2 = SynapseElife.getCamKII(nothing, res.xc, nothing)
    #append!(plasticity_obs.val.camk2,_camk2)
    _can = SynapseElife.getCaN(nothing, res.xc, nothing)
    #append!(plasticity_obs.val.can,_can)
    _weight = res.xd[38, :] .- res.xd[37, :]
    #append!(plasticity_obs.val.weight, _weight)

    plasticity_obs[] = (
        t = vcat(plasticity_obs[].t, res.time),
        camk2 = vcat(plasticity_obs[].camk2, _camk2),
        can = vcat(plasticity_obs[].can, _can),
        weight = vcat(plasticity_obs[].weight, _weight),
    )

    # update the progress bar
    next!(pbar; showvalues = [(:steps, length(events_sorted_times) + 1), (:t, p_synapse.t_end)])
    progress[] = 100
    notify(progress)

    @assert res.time[end] == p_synapse.t_end "Error in PDMP. Did not reach requested simulated time",
    "complete simulated time"

    return (t = tt, XC = XC, XD = XD)
end

export _evolve_synapse,
    firingPattern, dataProtocol, stp, initial_conditions_continuous_temp, initial_conditions_discrete
## blanket export of everything (less safe but quick for a utility file)
## At the end of utils_bonito.jl, before `end`:
# __all__ = [n for n in names(@__MODULE__; all = true) if !iskeyword(n) && n ∉ (:Main, :UtilsBonito, :Core, :Base)]
# export __all__...
end

# =============================================================================
#  SynapseElife Studio — Interactive GUI  (v6 — tabbed layout)
#  Built with Bonito (web UI) + WGLMakie (live plots) + CairoMakie (export)
#
#  ARCHITECTURE OVERVIEW
#  ─────────────────────
#  This file contains:
#    1. Helper functions   — input parsing, event construction
#    2. Simulation engine  — run_one_simulation! and run_simulation!
#    3. GUI helpers        — HTML table formatter
#    4. build_app()        — the full Bonito application (widgets, plots, callbacks)
#
#  TAB SYSTEM
#  ──────────
#  Bonito registers all widgets by their DOM id at page load time.
#  Hiding widgets with display:none prevents Bonito from finding them → crashes.
#  Solution: ALL widgets are rendered into the DOM once at startup (inside
#  invisible wrapper divs).  Tab switching is handled by a Julia-side
#  Observable{String} called `active_tab`.  The visible tab content is produced
#  by lift(active_tab) which returns a DOM node containing only the active
#  panel's content.  Inactive panels are rendered as empty divs.
#
#  SIMULATION MODES
#  ────────────────
#    • Protocol : reads parameters from a pre-defined experiment table
#    • STDP     : loops over Δt values, builds the full STDP curve
#    • HFS/LFS  : pre-synaptic train at a user-defined frequency
#    • TBS      : theta-burst stimulation (N bursts × P pulses @ intra Hz)
#    • Custom   : fully custom pre/post spike timestamps
#
#  EXPORT
#  ──────
#  WGLMakie renders in the browser via WebGL — it cannot write PNG files.
#  CairoMakie (server-side) is used exclusively for file export.
#  Tabular data and simulation results are exported as timestamped .xlsx files.
# =============================================================================

cd(@__DIR__)
using Pkg
pkg"activate ."

using Revise, SynapseElife,
Random,
WGLMakie,
Bonito,
PiecewiseDeterministicMarkovProcesses,
Sundials,
PrettyTables,
DataFrames,
CSV,
Statistics,
XLSX,
Dates

import CairoMakie

includet("utils_bonito.jl")
using .UtilsBonito


# =============================================================================
#  SECTION 1 — HELPER FUNCTIONS
# =============================================================================

# Parse a whitespace/comma-separated string of numbers into Vector{Float64}.
# Returns (values, error_message).  Empty string → empty vector (valid).
function parse_floats(s::String, name::String; allow_negative = true)
    cleaned = replace(strip(s), "," => " ")
    parts   = split(cleaned; keepempty = false)
    isempty(parts) && return Float64[], ""
    vals = Float64[]
    for p in parts
        v = tryparse(Float64, p)
        isnothing(v)               && return nothing, "⚠ $name: \"$p\" is not a valid number."
        (!allow_negative && v < 0) && return nothing, "⚠ $name: values must be ≥ 0."
        push!(vals, v)
    end
    return vals, ""
end

# Read and range-check a single Float64 from a TextField widget.
# Returns (value, error_message).
function get_valid_val(tb, lo, hi, name)
    v = tryparse(Float64, tb.value[])
    isnothing(v)       && return nothing, "⚠ $name: a number is required."
    (v < lo || v > hi) && return nothing, "⚠ $name must be between $lo and $hi."
    return v, ""
end

# Build the two vectors required by _evolve_synapse from absolute spike
# timestamps within one burst, repeated repeat_n times every repeat_interval ms.
# The output is chronologically sorted (mandatory for the PDMP integrator).
function build_events(pre_ts::Vector{Float64},
    post_ts::Vector{Float64},
    repeat_n::Int,
    repeat_interval::Float64,
    start_time::Float64)
    times  = Float64[]
    is_pre = Bool[]
    for rep in 0:(repeat_n - 1)
        offset = start_time + rep * repeat_interval
        for t in pre_ts;  push!(times, offset + t); push!(is_pre, true);  end
        for t in post_ts; push!(times, offset + t); push!(is_pre, false); end
    end
    perm = sortperm(times)
    return times[perm], is_pre[perm]
end

# ---------------------------------------------------------------------------
# estimate_memory_gb
# ---------------------------------------------------------------------------
# Estimates the RAM consumption of one _evolve_synapse call in GB.
# Based on the number of saved timepoints × number of continuous variables.
# Two contributions:
#   - GluON  : one short integration per pre-synaptic event (glu_width ms)
#   - GluOFF : one long integration per inter-event gap + final plateau
#
# n_pre      : number of pre-synaptic events
# t_end      : total simulation time (ms)
# glu_width  : glutamate pulse duration (ms), default 2.0 ms in SynapseElife
# n_xc       : number of continuous variables (default 40 — SynapseElife standard)
# sampling_rate : points per ms saved by the solver (default 0.01 → 1 pt/100ms)
function estimate_memory_gb(n_pre::Int, t_end::Float64;
                             glu_width     = 2.0,
                             n_xc          = 40,
                             sampling_rate = 0.01)
    # Points saved during GluON phases (one per pre event)
    pts_glu_on  = n_pre * glu_width * sampling_rate

    # Points saved during GluOFF phases (gaps between events + final plateau)
    # Conservative estimate: total time minus GluON time
    t_glu_off   = t_end - n_pre * glu_width
    pts_glu_off = max(0.0, t_glu_off) * sampling_rate

    total_pts   = pts_glu_on + pts_glu_off

    # Memory: each point stores n_xc Float64 values (8 bytes each)
    bytes = total_pts * n_xc * 8
    return bytes / 1e9  # convert to GB
end

# ---------------------------------------------------------------------------
# check_memory_safety
# ---------------------------------------------------------------------------
# Returns (ok, warning_message).
# ok = true  → safe to proceed
# ok = false → estimated RAM exceeds 2× available RAM, refuse to start
#
# Uses Sys.free_memory() for available RAM at call time.
# The 2× threshold accounts for Julia GC overhead and vector concatenation
# inside _evolve_synapse (vcat creates temporary copies).
function check_memory_safety(n_pre::Int, t_end::Float64; kwargs...)
    estimated_gb  = estimate_memory_gb(n_pre, t_end; kwargs...)
    available_gb  = Sys.free_memory() / 1e9
    total_gb      = Sys.total_memory() / 1e9
    threshold_gb  = 2 * available_gb

    msg = "Estimated RAM: $(round(estimated_gb; digits=2)) GB  |  " *
          "Available: $(round(available_gb; digits=2)) GB  |  " *
          "Total: $(round(total_gb; digits=2)) GB"
    println(msg)

    if estimated_gb > threshold_gb
        warn = "⚠ RAM warning: estimated $(round(estimated_gb; digits=1)) GB " *
               "> 2× available ($(round(threshold_gb; digits=1)) GB). " *
               "Reduce N pulses, N trains, or simulation duration."
        return false, warn
    end

    return true, ""
end

# =============================================================================
#  SECTION 2 — SIMULATION ENGINE
# =============================================================================

# Run one self-contained simulation and return the final synaptic weight (%).
# Used by the STDP loop: each Δt needs an independent run whose weight
# endpoint becomes one point on the STDP curve.
function run_one_simulation!(progress, plasticity_obs;
    pre_ts          = [0.0],
    post_ts         = Float64[],
    repeat_n        = 1,
    repeat_interval = 1000.0,   # ms between burst onsets
    ca_ext          = 2.5,      # mM  → converted to μM internally
    mg              = 1.0,
    temp            = 35.0,
    age             = "Adult")

    start   = 0.5e3
    _ca_ext = ca_ext * 1000.0
    _age    = (age == "Adult" ? 60.0 : 21.0)

    events_times, is_pre_or_post_event =
    build_events(pre_ts, post_ts, repeat_n, repeat_interval, start)

    afterend = 1.5e5   # 150 s margin after last event to observe the weight plateau
    t_end    = isempty(events_times) ? start + afterend : events_times[end] + afterend

    param_synapse = SynapseParams(
    t_end = t_end, soma_dist = 200.0, temp_rates = temp,
    Ca_ext = _ca_ext, Mg = mg, age = _age,
    injbap = 2.0, I_clamp = 1000.0, sampling_rate = 0.01)

    p           = param_synapse.p_release
    pre_synapse = PreSynapseParams(
    h = p[1] + p[2] / (1 + exp(p[3] * (_ca_ext - p[4]))))

    xc0 = initial_conditions_continuous_temp(param_synapse)
    xd0 = initial_conditions_discrete(param_synapse)

    # Short-term plasticity: which pre events actually release glutamate
    is_glu_release, _, _, _, _, bap_by_epsp_times =
    stp(param_synapse.t_end, pre_synapse, events_times, is_pre_or_post_event,
    _plot = false, algo = CHV(CVODE_BDF()))

    valid_idx = events_times .< param_synapse.t_end
    ev_clean  = events_times[valid_idx]
    ipp_clean = is_pre_or_post_event[valid_idx]

    # plasticity_obs is updated live inside _evolve_synapse so plots refresh
    # continuously during the simulation
    _evolve_synapse(
    xc0, xd0, param_synapse, ev_clean, ipp_clean,
    bap_by_epsp_times, is_glu_release,
    (CHV(CVODE_BDF(linear_solver = :GMRES)),
    CHV(CVODE_BDF(linear_solver = :GMRES)));
    abstol = 1e-6, reltol = 1e-5,
    save_positions = ((false, false), (false, true)),
    verbose = false, progress, plasticity_obs)

    return plasticity_obs[].weight[end]
end

# High-level wrapper: protocol mode reads from a DataFrame row;
# manual mode tiles pre_times/post_times at 1/freq Hz for `pulse` cycles.
function run_simulation!(progress, data, data_protocol, k;
    plasticity_obs,
    use_protocol = true,
    ca_ext = 2.5, mg = 1.0, temp = 35.0,
    pre_times  = [0.0],
    post_times = [10.0],
    pulse = 300, freq = 5.0, age = "Adult")
    _debug = true
    start  = 0.5e3

    if use_protocol
        @info "Protocol mode: $(data_protocol[k, :protocol])"
        pls           = data_protocol[!, :pulse][k]
        _freq         = data_protocol[!, :freq][k]
        _temp         = data_protocol[!, :temp][k]
        _ca_ext       = data_protocol[!, :exca][k]
        _mg           = data_protocol[!, :exmg][k]
        _age          = data_protocol[!, :age][k]
        _repeat_times = data_protocol[!, :repeat_times][k]
        _injbap       = data_protocol[!, :inj_time][k]
        _injection    = data_protocol[!, :injection][k]
        _AP_by_EPSP   = data_protocol[!, :AP_by_EPSP][k]

        events_times, is_pre_or_post_event = firingPattern(
        start_time   = start,
        n_pos        = data_protocol[!, :n_pos][k],
        delay_pos    = data_protocol[!, :delay_pos][k],
        n_pre        = data_protocol[!, :n_pre][k],
        delay_pre    = data_protocol[!, :delay_pre][k],
        delay        = data_protocol[!, :delay][k],
        pulse        = pls, freq = _freq,
        causal       = data_protocol[!, :causal][k],
        repeat_times = _repeat_times,
        repeat_after = data_protocol[!, :repeat_after][k])

        afterend = 1.5e5
        t_end    = start + 0.1 * (_repeat_times + 1) * pls * 1000 / _freq + afterend
    else
        @info "Manual / Custom pattern mode"
        pls    = pulse;  _freq = freq;  _temp = temp
        _ca_ext       = ca_ext * 1000.0
        _mg           = mg
        _age          = (age == "Adult" ? 60.0 : 21.0)
        _repeat_times = 0;  _injbap = 2.0;  _injection = 1000.0;  _AP_by_EPSP = "yes"

        # Tile the within-cycle pattern at the given frequency
        events_list = Tuple{Float64, Bool}[]
        for rep in 0:(pulse - 1)
            base_t = start + rep * (1000.0 / freq)
            for t_pre  in pre_times;  push!(events_list, (base_t + t_pre,  true));  end
            for t_post in post_times; push!(events_list, (base_t + t_post, false)); end
        end
        sort!(events_list; by = x -> x[1])
        events_times         = Float64[x[1] for x in events_list]
        is_pre_or_post_event = Bool[x[2]   for x in events_list]

        afterend = 1.5e5
        t_end    = start + 0.1 * pulse * 1000 / freq + afterend
    end

    param_synapse = SynapseParams(
    t_end = t_end, soma_dist = 200.0, temp_rates = _temp,
    Ca_ext = _ca_ext, Mg = _mg, age = _age,
    injbap = _injbap, I_clamp = _injection, sampling_rate = 0.01)

    p           = param_synapse.p_release
    pre_synapse = PreSynapseParams(
    h = p[1] + p[2] / (1 + exp(p[3] * (_ca_ext - p[4]))))

    xc0 = initial_conditions_continuous_temp(param_synapse)
    xd0 = initial_conditions_discrete(param_synapse)

    is_glu_release, _, _, _, _, bap_by_epsp_times =
    stp(param_synapse.t_end, pre_synapse, events_times, is_pre_or_post_event,
    _plot = false, algo = CHV(CVODE_BDF()))

    valid_idx = events_times .< param_synapse.t_end
    ev_clean  = events_times[valid_idx]
    ipp_clean = is_pre_or_post_event[valid_idx]
    _AP_by_EPSP_flag = use_protocol ? _AP_by_EPSP : "yes"

    _debug && @time _evolve_synapse(
    xc0, xd0, param_synapse, ev_clean, ipp_clean,
    _AP_by_EPSP_flag == "yes" ? bap_by_epsp_times : Float64[],
    is_glu_release,
    (CHV(CVODE_BDF(linear_solver = :GMRES)),
    CHV(CVODE_BDF(linear_solver = :GMRES)));
    abstol = 1e-6, reltol = 1e-5,
    save_positions = ((false, false), (false, true)),
    verbose = false, progress, plasticity_obs)

    use_protocol ? @warn(data_protocol[k, :protocol]) : @warn("Mode finished")
end


# =============================================================================
#  SECTION 3 — GUI HELPERS
# =============================================================================

# ---------------------------------------------------------------------------
# DOI mapping for Protocol tab
# ---------------------------------------------------------------------------
const ELIFE_DOI = "10.7554/eLife.80152"

const PROTOCOL_DOIS = Dict(
    "oconnor"                      => "10.1152/jn.00047",
    "Bittner"                      => "10.1126/science.aan3846",
    "Goldings01"                   => "10.1152/jn.2001.86.6.2998",
    "Buchenan"                     => "10.1113/jphysiol.2007.142984",
    "TigaretMellor16"              => "10.1038/ncomms10289",
    "TigaretMellor16_poisson"      => "10.1038/ncomms10289",
    "Tigaret_jitter"               => "10.1038/ncomms10289",
    "Tigaret_jitter_timespent"     => "10.1038/ncomms10289",
    "Tigaret_jitter_double_1"      => "10.1038/ncomms10289",
    "Tigaret_burst"                => "10.1038/ncomms10289",
    "Tigaret_burst_temp"           => "10.1038/ncomms10289",
    "Tigaret_burst_age"            => "10.1038/ncomms10289",
    "Tigaret_burst_freq"           => "10.1038/ncomms10289",
    "Tigaret_burst_ca"             => "10.1038/ncomms10289",
    "Tigaret_burst_dist"           => "10.1038/ncomms10289",
    "Tigaret_freq_1200"            => "10.1038/ncomms10289",
    "Tigaret_tripost"              => "10.1038/ncomms10289",
    "Tigaret_preburst"             => "10.1038/ncomms10289",
    "Tigaret_single"               => "10.1038/ncomms10289",
    "TigaretMellor_sparse"         => "10.1038/ncomms10289",
    "TigaretMellor_jitter_sparse"  => "10.1038/ncomms10289",
    "TigaretvsMeredith"            => "10.1038/ncomms10289",
    "DudekBear92-BCM"              => "10.1073/pnas.89.10.4363",
    "DudekBear92-BCM-900"          => "10.1073/pnas.89.10.4363",
    "DudekBear92-BCM-37"           => "10.1073/pnas.89.10.4363",
    "DudekBear92-BCM-33"           => "10.1073/pnas.89.10.4363",
    "DudekBear92-BCM-Ca"           => "10.1073/pnas.89.10.4363",
    "DudekBear92-BCM-priming"      => "10.1073/pnas.89.10.4363",
    "DudekBear92_timespend"        => "10.1073/pnas.89.10.4363",
    "DudekBear92-100"              => "10.1073/pnas.89.10.4363",
    "DudekBear_short"              => "10.1073/pnas.89.10.4363",
    "DudekBear92-sliding"          => "10.1073/pnas.89.10.4363",
    "DudekBear92_temp"             => "10.1073/pnas.89.10.4363",
    "DudekBear92_Ca"               => "10.1073/pnas.89.10.4363",
    "DudekBear92_Mg"               => "10.1073/pnas.89.10.4363",
    "DudekBear92_dist"             => "10.1073/pnas.89.10.4363",
    "DudekBear92-Age"              => "10.1073/pnas.89.10.4363",
    "DudekBear92_BCM_recovery"     => "10.1073/pnas.89.10.4363",
    "Dudek_jitter"                 => "10.1073/pnas.89.10.4363",
    "DudekBear93-LFS"              => "10.1523/JNEUROSCI.13-07-02910.1993",
    "DudekBear93-TBS"              => "10.1523/JNEUROSCI.13-07-02910.1993",
    "YannisDebanne20"              => "10.1113/jphysiol.1996.sp021204",
    "YannisDebanne20_freq"         => "10.1113/jphysiol.1996.sp021204",
    "YannisDebanne20_freq_delay"   => "10.1113/jphysiol.1996.sp021204",
    "YannisDebanne20_ratio"        => "10.1113/jphysiol.1996.sp021204",
    "YannisDebanne20_inv"          => "10.1113/jphysiol.1996.sp021204",
    "YannisDebanne_temp"           => "10.1113/jphysiol.1996.sp021204",
    "YannisDebanne_age"            => "10.1113/jphysiol.1996.sp021204",
    "YannisvsMeredith"             => "10.1113/jphysiol.1996.sp021204",
    "Cao-TBS"                      => "10.1152/jn.00780.2011",
    "Chang19"                      => "10.1038/s41467-019-10694-z",
    "Fujii_CaN"                    => "10.1016/j.celrep.2013.03.033",
    "Meredith03-GABA"              => "10.1523/JNEUROSCI.23-35-11142.2003",
    "WittenbergWang06_D"           => "10.1523/JNEUROSCI.5388-05.2006",
    "WittenbergWang06_B"           => "10.1523/JNEUROSCI.5388-05.2006",
    "WittenbergWang06_P"           => "10.1523/JNEUROSCI.5388-05.2006",
    "Mizuno01LTP"                  => "10.1046/j.0953-816x.2001.01679.x",
    "Mizuno01LTD"                  => "10.1046/j.0953-816x.2001.01679.x",
    "Mizuno01-LTP-Mg"              => "10.1046/j.0953-816x.2001.01679.x",
)

get_doi(protocol::String) = get(PROTOCOL_DOIS, protocol, ELIFE_DOI)

# Render a protocol DataFrame as a scrollable HTML table.
# The row at selected_idx is highlighted in blue.
function format_table(dp, selected_idx = 1)
    df = copy(dp)
    DataFrames.insertcols!(df, 1, "ID" => 1:DataFrames.nrow(df))

    if "paper" in names(df)
        df[!, :paper] = [get_doi(p) for p in df[!, :paper]]
        DataFrames.rename!(df, :paper => :DOI)
    end

    headers = [Bonito.DOM.th(string(name),
    style = "position:sticky;top:0;background-color:#2c3e50;" *
    "color:white;padding:8px;border:1px solid #ddd;z-index:10;")
    for name in names(df)]
    thead   = Bonito.DOM.thead(Bonito.DOM.tr(headers...))

    rows = map(enumerate(eachrow(df))) do (i, row)
        bg  = i == selected_idx ? "#dbeafe" : (iseven(i) ? "#f9f9f9" : "#ffffff")
        fw  = i == selected_idx ? "bold"    : "normal"
        col = i == selected_idx ? "#1e40af" : "black"
        cells = [Bonito.DOM.td(string(val),
        style = "padding:6px 8px;border:1px solid #ddd;white-space:nowrap;")
        for val in row]
        Bonito.DOM.tr(cells...;
        style = "background-color:$bg;font-weight:$fw;color:$col;transition:background-color 0.3s;")
    end
    return Bonito.DOM.table(thead, Bonito.DOM.tbody(rows...);
    style = "width:100%;border-collapse:collapse;font-family:'Segoe UI',sans-serif;font-size:12px;text-align:center;")
end

# ---------------------------------------------------------------------------
# run_batch!
# ---------------------------------------------------------------------------
# Runs the same simulation n_runs times independently and returns a vector
# of final weight change values (%). Each run gets fresh initial conditions
# and a new stochastic draw from the PDMP — this is the correct way to
# sample the distribution of plasticity outcomes.
#
# After each run, plasticity_obs[] is updated so the plots show the last
# run's dynamics. final_wc[] is updated with mean ± std after all runs.
# batch_results[] accumulates all weight values for CSV export.
#
# sim_fn: a zero-argument closure that runs one simulation and returns
#         the final weight. The closure captures all parameters from the
#         calling tab's Run callback.
function run_batch!(n_runs, p_bar, plasticity_obs, final_wc, batch_results,
                    all_run_curves, sim_fn)
    results = Float64[]

    all_run_curves[] = []

    for run in 1:n_runs
        plasticity_obs[] = (t=zeros(1), camk2=zeros(1), can=zeros(1), weight=zeros(1))
        p_bar[] = 0.0
        println("Run $run / $n_runs ...")

        w = sim_fn()
        push!(results, w)

        # Store curve for overlay (only when N ≤ 10 to avoid memory pressure)
        if n_runs <= 10
            obs = plasticity_obs[]
            all_run_curves[] = vcat(all_run_curves[],
                [(t = copy(obs.t), weight = copy(obs.weight))])
        end

        println("  Run $run => DW = $(round(w; digits=2))%")
    end

    batch_results[] = results

    if n_runs == 1
        final_wc[] = "Final weight change: $(round(results[1]; digits=1)) %"
    else
        m = mean(results); s = std(results)
        mn, mx = minimum(results), maximum(results)
        final_wc[] = "DW = $(round(m;digits=1)) ± $(round(s;digits=1)) %  " *
                     "[$(round(mn;digits=1)), $(round(mx;digits=1))]  N=$n_runs"
        println("Batch complete: mean=$(round(m;digits=2))% std=$(round(s;digits=2))%")
        if n_runs > 10
            println("N>10: showing last run curve. Mean DW in badge.")
        end
    end
end

# ---------------------------------------------------------------------------
# make_timeline_svg
# ---------------------------------------------------------------------------
# Generates an inline SVG timeline showing pre and post spike timestamps.
# Called inside tab_stdp_manual() whenever the Manual fields are shown.
# pre_times and post_times are vectors of Float64 timestamps in ms.
# The SVG is returned as a raw HTML string wrapped in a Bonito DOM node.
function make_timeline_svg(pre_times::Vector{Float64}, post_times::Vector{Float64})
    all_times = vcat(pre_times, post_times)
    isempty(all_times) && return Bonito.DOM.div()

    pre_str  = isempty(pre_times)  ? "—" : join(["$(round(Int,t)) ms" for t in pre_times],  "  ·  ")
    post_str = isempty(post_times) ? "—" : join(["$(round(Int,t)) ms" for t in post_times], "  ·  ")

    return Bonito.DOM.div(
        Bonito.DOM.p("Spike pattern preview (one cycle):";
            style = "font-size:11px;color:#7f8c8d;margin:10px 0 4px 0;"),
        Bonito.DOM.div(
            Bonito.DOM.span("▲ pre: ";  style="color:#2980b9;font-weight:bold;font-size:11px;"),
            Bonito.DOM.span(pre_str;    style="color:#2980b9;font-size:11px;"),
            Bonito.DOM.span("   |   ";  style="color:#bdc3c7;font-size:11px;"),
            Bonito.DOM.span("▼ post: "; style="color:#c0392b;font-weight:bold;font-size:11px;"),
            Bonito.DOM.span(post_str;   style="color:#c0392b;font-size:11px;");
            style = "background:#f8f9fa;border:1px solid #ecf0f1;border-radius:6px;" *
                    "padding:8px 12px;display:flex;align-items:center;flex-wrap:wrap;gap:4px;");
        style = "margin-top:12px;")
end

# =============================================================================
#  SECTION 4 — BONITO APPLICATION
# =============================================================================

function build_app(; fontsize = 20)
    app = Bonito.App() do session

        # ── SHARED OBSERVABLES ────────────────────────────────────────────────
        p_bar          = Observable{Float64}(0.0)
        plasticity_obs = Observable((t=zeros(1), camk2=zeros(1),
        can=zeros(1), weight=zeros(1)))
        data             = Observable(Float64[])
        error_msg        = Observable("")
        confirm_stop_msg = Observable("")
        current_task     = Ref{Task}()
        final_wc = Observable("—") # stores the final weight change as a formatted string
        batch_results = Observable(Float64[])  # stores all weight values from batch runs
        stdp_points   = Observable((dt = Float64[], dw = Float64[]))
        all_run_curves = Observable(Vector{NamedTuple{(:t, :weight),
                            Tuple{Vector{Float64}, Vector{Float64}}}}())


        # active_tab: the single source of truth for which tab is shown.
        # Possible values:| "hfs" | "tbs" | "stdp/manual" | "protocol"
        # Default is "stdp".
        active_tab = Observable("hfs")

        # ── GLOBAL TOOLBAR WIDGETS ────────────────────────────────────────────
        btn_run  = Bonito.Button("▶  Run")
        btn_stop = Bonito.Button("■  Stop")
        btn_era  = Bonito.Button("✕  Erase")
        btn_csv  = Bonito.Button("Export Excel")
        btn_plt  = Bonito.Button("Export Plots")
        tb_n_runs = Bonito.TextField("1")  # number of batch runs (default 1 = current behaviour)


        # ── TAB BUTTONS ───────────────────────────────────────────────────────
        # Each button sets active_tab[] to the corresponding tab name.
        # They are plain Bonito.Button widgets — no JS onclick needed.
        btn_tab_hfs      = Bonito.Button("HFS / LFS")
        btn_tab_tbs      = Bonito.Button("TBS")
        btn_tab_stdp_manual = Bonito.Button("Custom Stimulation")
        btn_tab_protocol = Bonito.Button("Protocol")

        # ── PROTOCOL TAB WIDGETS ──────────────────────────────────────────────
        protocol_options  = collect(SynapseElife._protocols)
        idx_default       = something(findfirst(==("TigaretMellor16"), protocol_options), 1)
        dp_init           = dataProtocol(protocol_options[idx_default])
        menu_prot         = Bonito.Dropdown(protocol_options; index = idx_default)
        menu_exp          = Bonito.Dropdown(collect(1:size(dp_init, 1)))
        data_protocol_str = Observable(format_table(dp_init, 1))

        # ── MANUAL TAB WIDGETS ────────────────────────────────────────────────
        # Pre times / Post times: space-separated absolute timestamps (ms)
        # within one burst cycle.  Example: Pre "0 10 20", Post "30".
        tb_man_pre    = Bonito.TextField("0")
        tb_man_post   = Bonito.TextField("10")
        tb_man_pls    = Bonito.TextField("300")   # number of cycle repetitions
        tb_man_freq   = Bonito.TextField("5.0")   # repetition frequency (Hz)
        tb_man_ca     = Bonito.TextField("2.5")
        tb_man_mg     = Bonito.TextField("1.0")
        tb_man_temp   = Bonito.TextField("35.0")
        menu_man_age  = Bonito.Dropdown(["Young", "Adult"]; index = 2)

        # ── HFS / LFS TAB WIDGETS ─────────────────────────────────────────────
        # Pure pre-synaptic train (no post spikes).
        # Label HFS vs LFS is determined automatically: ≥ 50 Hz → HFS.
        tb_hfs_freq   = Bonito.TextField("100.0")
        tb_hfs_pulses = Bonito.TextField("100")
        tb_hfs_ntr    = Bonito.TextField("1")
        tb_hfs_itr    = Bonito.TextField("15000.0")
        tb_hfs_ca     = Bonito.TextField("2.5")
        tb_hfs_mg     = Bonito.TextField("1.0")
        tb_hfs_temp   = Bonito.TextField("35.0")
        menu_hfs_age  = Bonito.Dropdown(["Young", "Adult"]; index = 2)

        # ── TBS TAB WIDGETS ───────────────────────────────────────────────────
        # Classic TBS protocol: N bursts × P pulses/burst at freq_intra Hz,
        # repeated at freq_theta Hz (theta period = 1000/freq_theta ms).
        # Intra-burst spike offsets: [0, dt, 2dt, ...] where dt = 1000/freq_intra ms.
        tb_tbs_nb     = Bonito.TextField("10")     # number of bursts per train
        tb_tbs_ppb    = Bonito.TextField("4")      # pulses per burst
        tb_tbs_fi     = Bonito.TextField("100.0")  # intra-burst frequency (Hz)
        tb_tbs_ft     = Bonito.TextField("5.0")    # theta / inter-burst frequency (Hz)
        tb_tbs_ntr    = Bonito.TextField("10")     # number of trains
        tb_tbs_itr    = Bonito.TextField("15000.0") # inter-train interval (ms)
        tb_tbs_ca     = Bonito.TextField("2.5")
        tb_tbs_mg     = Bonito.TextField("1.0")
        tb_tbs_temp   = Bonito.TextField("35.0")
        menu_tbs_age  = Bonito.Dropdown(["Young", "Adult"]; index = 2)

        # ── WGLMAKIE FIGURES ──────────────────────────────────────────────────

        # Main figure: 3 plots + progress bar
        fig = Figure(; size = (1100, 650), fontsize)

        ax_camk = Axis(fig[1, 1]; xlabel = "time (s)", ylabel = "μM",
        title = "CaMKII & CaN",
        ylabelsize = fontsize/1.5, xlabelsize = fontsize/1.5,
        yticklabelsize = fontsize/2, xticklabelsize = fontsize/2)
        ylims!(ax_camk, 0, 34); xlims!(ax_camk, 0, 200)

        ax_phase = Axis(fig[1, 2]; ylabel = "CaMKII (μM)", xlabel = "CaN (μM)",
        title = "Phase portrait",
        ylabelsize = fontsize/1.5, xlabelsize = fontsize/1.5,
        yticklabelsize = fontsize/2, xticklabelsize = fontsize/2)

        ax_wt = Axis(fig[1, 3]; ylabel = "weight change (%)", xlabel = "time (s)",
        title = "Synaptic weight",
        ylabelsize = fontsize/1.5, xlabelsize = fontsize/1.5,
        yticklabelsize = fontsize/2, xticklabelsize = fontsize/2)
        xlims!(ax_wt, 0, 200); ylims!(ax_wt, -150, 150)
        wt_main_color    = :steelblue
        wt_overlay_alpha = 0.45

        ax_prog = Axis(fig[2, 1:3]; title = "Progress", height = 16, yticks = [0])
        pbar_rect = map(p_bar) do x; Rect(0.0, 0.0, x, 1.0); end
        poly!(ax_prog, Rect(0.0, 0.0, 100.0, 1.0);
        color = :white, strokecolor = :red, strokewidth = 1)
        poly!(ax_prog, pbar_rect; color = :red)

        # Live curves (lift makes them reactive to plasticity_obs)
        lines!(ax_camk, lift(x -> x.t./1000, plasticity_obs),
        lift(x -> x.camk2,   plasticity_obs); label = "CaMKII")
        lines!(ax_camk, lift(x -> x.t./1000, plasticity_obs),
        lift(x -> x.can,     plasticity_obs); label = "CaN")
        axislegend(ax_camk)

        lines!(ax_phase, lift(x -> x.can,   plasticity_obs),
        lift(x -> x.camk2, plasticity_obs))

        # Current run — always shown in steelblue
        lines!(ax_wt, lift(x -> x.t./1000, plasticity_obs),
               lift(x -> x.weight, plasticity_obs);
               color = wt_main_color, linewidth = 2)

        # Overlay: each completed run adds one thin line when all_run_curves updates
        on(all_run_curves) do curves
            n = length(curves)
            n < 2 && return
            c = curves[end]
            α = Float32((n - 1) / max(9, n - 1))   # 0 → 1 across up to 10 runs
            col = RGBAf(0.16f0 + 0.04f0*α, 0.50f0 - 0.15f0*α, 0.73f0 - 0.20f0*α,
                        wt_overlay_alpha)
            lines!(ax_wt, c.t ./ 1000, c.weight; color = col, linewidth = 1.0f0)
        end

        poly!(ax_phase, Point2f[(P[1],P[2]) for P in SynapseParams().LTP_region.vertices];
        color = ("red",  0.1))
        poly!(ax_phase, Point2f[(P[1],P[2]) for P in SynapseParams().LTD_region.vertices];
        color = ("blue", 0.1))


        # ── CSS CONSTANTS ─────────────────────────────────────────────────────
        # Reusable style strings to keep the DOM section readable
        CARD   = "background:white;padding:22px;border-radius:10px;" *
        "box-shadow:0 3px 10px rgba(0,0,0,0.09);margin-bottom:0;"
        TITLE  = "margin-top:0;margin-bottom:14px;color:#2c3e50;" *
        "font-family:'Segoe UI',sans-serif;font-size:17px;" *
        "border-bottom:2px solid #ecf0f1;padding-bottom:7px;"
        LBL    = "font-size:12px;color:#7f8c8d;"
        COL    = "display:flex;flex-direction:column;gap:4px;"
        GRID   = "display:grid;grid-template-columns:repeat(auto-fit,minmax(145px,1fr));gap:14px;"
        HINT   = "font-size:11px;color:#2c3e50;margin-bottom:12px;line-height:1.5;"
        # Label + circular ? badge with tooltip on hover.
        function lbl_help(text::String; tooltip::String="", bold=false)
            fw = bold ? "font-weight:bold;" : ""
            Bonito.DOM.div(
                Bonito.DOM.span(text; style=LBL * fw),
                isempty(tooltip) ? Bonito.DOM.span() :
                Bonito.DOM.span("?";
                    style = "display:inline-flex;align-items:center;justify-content:center;" *
                            "width:14px;height:14px;border-radius:50%;background:#bdc3c7;" *
                            "color:white;font-size:9px;font-weight:bold;margin-left:4px;" *
                            "cursor:help;flex-shrink:0;",
                    title = tooltip);
                style = "display:flex;align-items:center;gap:2px;")
        end
        # Positional compat: lbl_help("label", "tooltip"; bold=false)
        lbl_help(text::String, tooltip::String; bold=false) =
            lbl_help(text; tooltip=tooltip, bold=bold)

        # Helper: 4 shared physical-parameter fields used by all stimulation tabs
        function phys_fields(tb_c, tb_m, tb_t, menu_a)
            Bonito.DOM.div(
            Bonito.DOM.div(lbl_help("Ca²⁺ [2.0–4] (mM):",
                "Extracellular calcium concentration (mM). Standard: 2.5 mM. Higher values increase release probability and LTP."), tb_c, style=COL),
            Bonito.DOM.div(lbl_help("Mg²⁺ [0.8–3] (mM):",
                "Extracellular magnesium concentration (mM). Mg²⁺ blocks NMDA receptors at rest. Higher Mg²⁺ reduces NMDA activation."), tb_m, style=COL),
            Bonito.DOM.div(lbl_help("Temp [20–40] (°C):",
                "Bath temperature (°C). Affects reaction rates and ion channel kinetics. Typical in vitro: 32–35°C."), tb_t, style=COL),
            Bonito.DOM.div(lbl_help("Age:",
                "Animal age. Young = 21 days (juvenile rat), Adult = 60 days. Affects release probability and synaptic parameters."), menu_a, style=COL);
            style = GRID)
        end

        # ── TAB CONTENT BUILDERS ──────────────────────────────────────────────
        # Each function returns a DOM node for that tab's parameter panel.
        # They are called inside lift(active_tab) so only the active one
        # is rendered in the visible part of the page.

        function tab_stdp_manual()

            # ── Live timeline preview ──────────────────────────────────────────────
            timeline = lift(tb_man_pre.value, tb_man_post.value) do pre_str, post_str
                pre_vals,  _ = parse_floats(pre_str,  "pre";  allow_negative=false)
                post_vals, _ = parse_floats(post_str, "post"; allow_negative=false)
                pre_vals  = isnothing(pre_vals)  ? Float64[] : pre_vals
                post_vals = isnothing(post_vals) ? Float64[] : post_vals
                make_timeline_svg(pre_vals, post_vals)
            end

            Bonito.DOM.div(
                Bonito.DOM.h3("Custom Stimulation"; style = TITLE),

                # ── Description + exemples ─────────────────────────────────────────
                Bonito.DOM.div(
                    Bonito.DOM.p(
                        "Define any pre/post-synaptic spike pattern by entering absolute timestamps (ms) " *
                        "within one cycle. The pattern is then repeated N times at the given frequency.";
                        style = HINT),

                    Bonito.DOM.div(
                        Bonito.DOM.p("STDP examples:";
                            style = "font-size:11px;font-weight:bold;color:#2c3e50;margin:0 0 4px 0;"),
                        Bonito.DOM.div(
                            Bonito.DOM.span("▲ Pre → Post (LTP expected): ";
                                style="font-size:11px;color:#2980b9;font-weight:bold;"),
                            Bonito.DOM.span("Pre = 0  |  Post = 10";
                                style="font-size:11px;color:#2c3e50;font-family:monospace;"),
                            Bonito.DOM.span("  → Δt = +10 ms";
                                style="font-size:11px;color:#7f8c8d;");
                            style="display:flex;align-items:center;gap:6px;margin-bottom:3px;"),
                        Bonito.DOM.div(
                            Bonito.DOM.span("▼ Post → Pre (LTD expected): ";
                                style="font-size:11px;color:#c0392b;font-weight:bold;"),
                            Bonito.DOM.span("Pre = 10  |  Post = 0";
                                style="font-size:11px;color:#2c3e50;font-family:monospace;"),
                            Bonito.DOM.span("  → Δt = −10 ms";
                                style="font-size:11px;color:#7f8c8d;");
                            style="display:flex;align-items:center;gap:6px;margin-bottom:3px;"),
                        Bonito.DOM.div(
                            Bonito.DOM.span(" Multi-spike: ";
                                style="font-size:11px;color:#8e44ad;font-weight:bold;"),
                            Bonito.DOM.span("Pre = 0 10 20  |  Post = 5 15";
                                style="font-size:11px;color:#2c3e50;font-family:monospace;"),
                            Bonito.DOM.span("  → 3 pre + 2 post per cycle";
                                style="font-size:11px;color:#7f8c8d;");
                            style="display:flex;align-items:center;gap:6px;margin-bottom:3px;"),
                        Bonito.DOM.div(
                            Bonito.DOM.span("● Pure pre-synaptic (HFS/LFS): ";
                                style="font-size:11px;color:#1a6b3c;font-weight:bold;"),
                            Bonito.DOM.span("Pre = 0  |  Post = (empty)";
                                style="font-size:11px;color:#2c3e50;font-family:monospace;");
                            style="display:flex;align-items:center;gap:6px;");
                        style = "background:#f8f9fa;border-radius:6px;padding:10px 14px;" *
                                "border-left:3px solid #bdc3c7;margin-bottom:14px;");

                    style = "margin-bottom:4px;"),

                # ── Champs de saisie ───────────────────────────────────────────────
                Bonito.DOM.div(
                    Bonito.DOM.p("Stimulation";
                        style="font-size:11px;font-weight:bold;color:#2c3e50;margin:0 0 8px 0;"),
                    Bonito.DOM.div(
                        Bonito.DOM.div(lbl_help("Pre times (ms):";
                            tooltip = "Absolute timestamps of pre-synaptic spikes within one cycle (ms), " *
                                      "separated by spaces or commas. " *
                                      "Example: '0 10 20' → 3 spikes at t=0, t=10, t=20 ms. " *
                                      "These timestamps are repeated N times at the given frequency.",
                            bold=true), tb_man_pre, style=COL),
                        Bonito.DOM.div(lbl_help("Post times (ms):";
                            tooltip = "Absolute timestamps of post-synaptic spikes within one cycle (ms), " *
                                      "separated by spaces or commas. " *
                                      "Example: '5' → 1 spike at t=5 ms. " *
                                      "Leave empty for pure pre-synaptic stimulation (HFS/LFS-style).",
                            bold=true), tb_man_post, style=COL),
                        Bonito.DOM.div(lbl_help("Repetitions:";
                            tooltip = "Number of times the pre/post pattern is repeated. " *
                                      "Each repetition starts at t = rep × (1000 / frequency) ms. " *
                                      "For STDP: use 300 repetitions at 5 Hz (Tigaret et al. 2016)."),
                            tb_man_pls, style=COL),
                        Bonito.DOM.div(lbl_help("Frequency (Hz):";
                            tooltip = "Repetition frequency of the pattern (Hz). " *
                                      "5 Hz → 200 ms between cycles. 10 Hz → 100 ms between cycles."),
                            tb_man_freq, style=COL);
                        style = GRID);
                    style = "background:#f8f9fa;border-radius:8px;padding:14px;margin-bottom:10px;"),

                Bonito.DOM.div(
                    Bonito.DOM.p("Main parameters";
                        style="font-size:11px;font-weight:bold;color:#2c3e50;margin:0 0 8px 0;"),
                    phys_fields(tb_man_ca, tb_man_mg, tb_man_temp, menu_man_age);
                    style = "background:#f8f9fa;border-radius:8px;padding:14px;margin-bottom:10px;"),

                # ── Timeline preview ───────────────────────────────────────────────
                Bonito.DOM.div(timeline);

                style = CARD * "border-left:4px solid #6c3483;")
        end

        function tab_hfs()
            Bonito.DOM.div(
                Bonito.DOM.h3("High / Low Frequency Stimulation"; style=TITLE),
                Bonito.DOM.p(
                    "Simulates a pure pre-synaptic spike train (no post-synaptic spikes). " *
                    "Frequency ≥ 50 Hz → HFS; < 50 Hz → LFS. " *
                    "Dudek & Bear 1992: 900 pulses at 1/3/5/10/50 Hz, 1 train.";
                    style = HINT),

                # ── Boîte Stimulation ──────────────────────────────────────────────
                Bonito.DOM.div(
                    Bonito.DOM.p("Stimulation";
                        style="font-size:11px;font-weight:bold;color:#2c3e50;margin:0 0 8px 0;"),
                    Bonito.DOM.div(
                        Bonito.DOM.div(lbl_help("Frequency (Hz):",
                            "Stimulation frequency in Hz. ≥ 50 Hz → HFS (induces LTP). < 50 Hz → LFS (may induce LTD). Dudek & Bear 1992: 1, 3, 5, 10, 50 Hz."),
                            tb_hfs_freq, style=COL),
                        Bonito.DOM.div(lbl_help("Pulses / train:",
                            "Number of pre-synaptic spikes per train. Dudek & Bear 1992: 900 pulses per train."),
                            tb_hfs_pulses, style=COL),
                        Bonito.DOM.div(lbl_help("N trains:",
                            "Number of times the full train is repeated. N=1: single continuous train. N>1: trains separated by the inter-train interval."),
                            tb_hfs_ntr, style=COL),
                        Bonito.DOM.div(lbl_help("Inter-train (ms):",
                            "Pause between trains in milliseconds. Only used when N trains > 1. Example: 15000 ms = 15 s."),
                            tb_hfs_itr, style=COL);
                        style = GRID);
                    style = "background:#f8f9fa;border-radius:8px;padding:14px;margin-bottom:10px;"),

                # ── Boîte Main parameters ───────────────────────────────────────────
                Bonito.DOM.div(
                    Bonito.DOM.p("Main parameters";
                        style="font-size:11px;font-weight:bold;color:#2c3e50;margin:0 0 8px 0;"),
                    phys_fields(tb_hfs_ca, tb_hfs_mg, tb_hfs_temp, menu_hfs_age);
                    style = "background:#f8f9fa;border-radius:8px;padding:14px;");

                style = CARD * "border-left:4px solid #1a6b3c;")
        end

        function tab_tbs()
            Bonito.DOM.div(
                Bonito.DOM.h3("Theta Burst Stimulation"; style=TITLE),
                Bonito.DOM.p(
                    "Classic TBS: N trains × N bursts × P pulses/burst at intra-burst frequency, " *
                    "repeated at theta frequency, with a pause between trains. " *
                    "Houeland et al. 2010: 10 trains × 10 bursts × 4 pulses @ 100 Hz, " *
                    "theta = 5 Hz, inter-train = 15 s.";
                    style = HINT),

                # ── Boîte Stimulation ──────────────────────────────────────────────
                Bonito.DOM.div(
                    Bonito.DOM.p("Stimulation";
                        style="font-size:11px;font-weight:bold;color:#2c3e50;margin:0 0 8px 0;"),
                    Bonito.DOM.div(
                        Bonito.DOM.div(lbl_help("N trains:",
                            "Number of TBS trains. Houeland et al. 2010: 10 trains."), tb_tbs_ntr, style=COL),
                        Bonito.DOM.div(lbl_help("Inter-train (ms):",
                            "Pause between trains in milliseconds. Houeland et al. 2010: 15000 ms (15 s)."), tb_tbs_itr, style=COL),
                        Bonito.DOM.div(lbl_help("N bursts / train:",
                            "Number of bursts per train. Each burst contains several pulses at high frequency. Houeland et al. 2010: 10 bursts per train."), tb_tbs_nb, style=COL),
                        Bonito.DOM.div(lbl_help("Pulses / burst:",
                            "Number of pre-synaptic spikes within one burst. Houeland et al. 2010: 4 pulses per burst."), tb_tbs_ppb, style=COL),
                        Bonito.DOM.div(lbl_help("Intra freq (Hz):",
                            "Frequency of spikes within one burst. Houeland et al. 2010: 100 Hz (10 ms between spikes)."), tb_tbs_fi, style=COL),
                        Bonito.DOM.div(lbl_help("Theta freq (Hz):",
                            "Frequency between burst onsets within one train. Typical theta rhythm: 5 Hz (200 ms between bursts)."), tb_tbs_ft, style=COL);
                        style = GRID);
                    style = "background:#f8f9fa;border-radius:8px;padding:14px;margin-bottom:10px;"),

                # ── Boîte Main parameters ───────────────────────────────────────────
                Bonito.DOM.div(
                    Bonito.DOM.p("Main parameters";
                        style="font-size:11px;font-weight:bold;color:#2c3e50;margin:0 0 8px 0;"),
                    phys_fields(tb_tbs_ca, tb_tbs_mg, tb_tbs_temp, menu_tbs_age);
                    style = "background:#f8f9fa;border-radius:8px;padding:14px;");

                style = CARD * "border-left:4px solid #b45309;")
        end

        function tab_protocol()
            Bonito.DOM.div(
            Bonito.DOM.h3("Protocol Mode"; style=TITLE),
            Bonito.DOM.p(
            "Protocols from Rodrigues et al. 2023 (eLife, DOI: 10.7554/eLife.80152). " *
            "Select a protocol and an experiment. " *
            "All simulation parameters are read from the table row. " *
            "The selected row is highlighted in blue.";
            style = HINT),
            Bonito.DOM.div(
            Bonito.DOM.span("Protocol:", style="font-weight:bold;font-size:13px;"), menu_prot,
            Bonito.DOM.span("  Exp #:",  style="font-weight:bold;font-size:13px;margin-left:10px;"), menu_exp;
            style="display:flex;gap:10px;align-items:center;margin-bottom:14px;"),
            Bonito.DOM.div(data_protocol_str;
            style="overflow-y:auto;overflow-x:auto;max-height:260px;" *
            "border:1px solid #ecf0f1;border-radius:8px;");
            style = CARD * "border-left:4px solid #2980b9;")
        end

        # ── CALLBACKS ─────────────────────────────────────────────────────────

        # Tab button callbacks — each sets active_tab to the tab's key string
        on(session, btn_tab_hfs.value)      do _; active_tab[] = "hfs";      end
        on(session, btn_tab_tbs.value)      do _; active_tab[] = "tbs";      end
        on(session, btn_tab_stdp_manual.value) do _; active_tab[] = "stdp_manual"; end
        on(session, btn_tab_protocol.value) do _; active_tab[] = "protocol"; end

        # Protocol table updates
        on(session, menu_prot.value) do s
            dp = dataProtocol(s)
            menu_exp.options[]  = collect(1:size(dp, 1))
            menu_exp.value[]    = 1
            data_protocol_str[] = format_table(dp, 1)
        end
        on(session, menu_exp.value) do idx
            isnothing(idx) && return
            data_protocol_str[] = format_table(dataProtocol(menu_prot.value[]), idx)
        end

        # Erase: reset plots and progress bar (does not stop a running simulation)
        on(session, btn_era.value) do _
            p_bar[]          = 0.0
            plasticity_obs[] = (t=zeros(1), camk2=zeros(1),
            can=zeros(1), weight=zeros(1))
            final_wc[]       = "—"
            all_run_curves[] = []
            empty!(ax_wt)
            lines!(ax_wt, lift(x -> x.t./1000, plasticity_obs),
                          lift(x -> x.weight,  plasticity_obs);
                   color = wt_main_color, linewidth = 2)
            xlims!(ax_wt, 0, 200); ylims!(ax_wt, -150, 150)
        end

        # Stop: interrupt the running async task; plots are preserved
        on(session, btn_stop.value) do _
            if isassigned(current_task) && !istaskdone(current_task[])
                Base.throwto(current_task[], InterruptException())
                confirm_stop_msg[] = "Simulation stopped — plots preserved."
            else
                confirm_stop_msg[] = "No simulation running."
            end
            @async (sleep(3); confirm_stop_msg[] = "")
        end

        # Export Excel (.xlsx) — timestamped, two sheets: Results + Configuration
        on(session, btn_csv.value) do _
            obs       = plasticity_obs[]
            timestamp = Dates.format(Dates.now(), "yyyy-mm-dd_HH-MM-SS")
            filename  = "SynapseElife_$(timestamp).xlsx"

            XLSX.openxlsx(filename, mode="w") do xf

                # ── Feuille 1 : Results ───────────────────────────────────────────
                sheet1 = xf[1]
                XLSX.rename!(sheet1, "Results")

                if !isempty(stdp_points[].dt)
                    # Mode STDP : une ligne par Δt
                    sheet1["A1"] = "Mode"
                    sheet1["B1"] = "STDP"
                    sheet1["A3"] = "Delta_t_ms"
                    sheet1["B3"] = "WeightChange_pct"
                    for (i, (dt, dw)) in enumerate(zip(stdp_points[].dt, stdp_points[].dw))
                        sheet1["A$(i+3)"] = dt
                        sheet1["B$(i+3)"] = dw
                    end
                    println("XLSX → STDP results ($(length(stdp_points[].dt)) points)")

                elseif !isempty(batch_results[])
                    # Mode Batch : une ligne par run
                    sheet1["A1"] = "Mode"
                    sheet1["B1"] = "Batch"
                    sheet1["A3"] = "Run"
                    sheet1["B3"] = "WeightChange_pct"
                    for (i, w) in enumerate(batch_results[])
                        sheet1["A$(i+3)"] = i
                        sheet1["B$(i+3)"] = w
                    end
                    # Stats en bas
                    n = length(batch_results[])
                    m = mean(batch_results[]); s = std(batch_results[])
                    sheet1["A$(n+5)"] = "Mean"
                    sheet1["B$(n+5)"] = m
                    sheet1["A$(n+6)"] = "Std"
                    sheet1["B$(n+6)"] = s
                    sheet1["A$(n+7)"] = "Min"
                    sheet1["B$(n+7)"] = minimum(batch_results[])
                    sheet1["A$(n+8)"] = "Max"
                    sheet1["B$(n+8)"] = maximum(batch_results[])
                    println("XLSX → Batch results ($n runs)")

                else
                    # Mode single run : série temporelle complète
                    sheet1["A1"] = "Mode"
                    sheet1["B1"] = "Single run"
                    sheet1["A3"] = "Time_s"
                    sheet1["B3"] = "CaMKII_uM"
                    sheet1["C3"] = "CaN_uM"
                    sheet1["D3"] = "WeightChange_pct"
                    t_sec = obs.t ./ 1000
                    for i in 1:length(t_sec)
                        sheet1["A$(i+3)"] = t_sec[i]
                        sheet1["B$(i+3)"] = obs.camk2[i]
                        sheet1["C$(i+3)"] = obs.can[i]
                        sheet1["D$(i+3)"] = obs.weight[i]
                    end
                    println("XLSX → Single run timeseries ($(length(t_sec)) points)")
                end

                # ── Feuille 2 : Configuration ─────────────────────────────────────
                XLSX.addsheet!(xf, "Configuration")
                sheet2 = xf[2]

                sheet2["A1"] = "Parameter"
                sheet2["B1"] = "Value"

                # Paramètres communs
                params = [
                    ("Timestamp",     timestamp),
                    ("Active tab",    active_tab[]),
                    ("N runs",        tb_n_runs.value[]),
                    ("Final WC",      final_wc[]),
                ]

                # Paramètres spécifiques à l'onglet actif
                tab = active_tab[]
                if tab == "hfs"
                    append!(params, [
                        ("Frequency Hz",   tb_hfs_freq.value[]),
                        ("Pulses / train", tb_hfs_pulses.value[]),
                        ("N trains",       tb_hfs_ntr.value[]),
                        ("Inter-train ms", tb_hfs_itr.value[]),
                        ("Ca2+ mM",        tb_hfs_ca.value[]),
                        ("Mg2+ mM",        tb_hfs_mg.value[]),
                        ("Temp C",         tb_hfs_temp.value[]),
                        ("Age",            menu_hfs_age.value[]),
                    ])
                elseif tab == "tbs"
                    append!(params, [
                        ("N trains",       tb_tbs_ntr.value[]),
                        ("Inter-train ms", tb_tbs_itr.value[]),
                        ("N bursts/train", tb_tbs_nb.value[]),
                        ("Pulses/burst",   tb_tbs_ppb.value[]),
                        ("Intra freq Hz",  tb_tbs_fi.value[]),
                        ("Theta freq Hz",  tb_tbs_ft.value[]),
                        ("Ca2+ mM",        tb_tbs_ca.value[]),
                        ("Mg2+ mM",        tb_tbs_mg.value[]),
                        ("Temp C",         tb_tbs_temp.value[]),
                        ("Age",            menu_tbs_age.value[]),
                    ])
                elseif tab == "stdp_manual"
                    append!(params, [
                        ("Pre times ms",   tb_man_pre.value[]),
                        ("Post times ms",  tb_man_post.value[]),
                        ("Repetitions",    tb_man_pls.value[]),
                        ("Frequency Hz",   tb_man_freq.value[]),
                        ("Ca2+ mM",        tb_man_ca.value[]),
                        ("Mg2+ mM",        tb_man_mg.value[]),
                        ("Temp C",         tb_man_temp.value[]),
                        ("Age",            menu_man_age.value[]),
                    ])
                elseif tab == "protocol"
                    k_str = string(menu_exp.value[])
                    append!(params, [
                        ("Protocol",       menu_prot.value[]),
                        ("Experiment #",   k_str),
                    ])
                end

                for (i, (k, v)) in enumerate(params)
                    sheet2["A$(i+1)"] = string(k)
                    sheet2["B$(i+1)"] = string(v)
                end

            end  # close XLSX file

            println("Exported → $filename")
        end

        # Export Plots (CairoMakie — server-side PNG rendering)
        on(session, btn_plt.value) do _
            obs   = plasticity_obs[]
            t_sec = obs.t ./ 1000
            ltp   = SynapseParams().LTP_region.vertices
            ltd   = SynapseParams().LTD_region.vertices

            f1 = CairoMakie.Figure(size=(800,500))
            a1 = CairoMakie.Axis(f1[1,1]; xlabel="Time (s)", ylabel="Concentration (μM)")
            CairoMakie.lines!(a1, t_sec, obs.camk2; label="CaMKII", linewidth=2)
            CairoMakie.lines!(a1, t_sec, obs.can;   label="CaN",    linewidth=2)
            CairoMakie.axislegend(a1)
            CairoMakie.save("plot_1_enzymes.png", f1; px_per_unit=2)

            f2 = CairoMakie.Figure(size=(600,500))
            a2 = CairoMakie.Axis(f2[1,1]; xlabel="CaN (μM)", ylabel="CaMKII (μM)")
            CairoMakie.poly!(a2, CairoMakie.Point2f[(P[1],P[2]) for P in ltp]; color=(:red,0.1))
            CairoMakie.poly!(a2, CairoMakie.Point2f[(P[1],P[2]) for P in ltd]; color=(:blue,0.1))
            CairoMakie.lines!(a2, obs.can, obs.camk2; linewidth=2)
            CairoMakie.save("plot_2_phase.png", f2; px_per_unit=2)

            f3 = CairoMakie.Figure(size=(800,500))
            a3 = CairoMakie.Axis(f3[1,1]; xlabel="Time (s)", ylabel="Weight change (%)")
            CairoMakie.lines!(a3, t_sec, obs.weight; linewidth=2)
            CairoMakie.save("plot_3_weight.png", f3; px_per_unit=2)

            if !isempty(stdp_points[].dt)
                sp = stdp_points[]
                f4 = CairoMakie.Figure(size=(600,450))
                a4 = CairoMakie.Axis(f4[1,1]; xlabel="Δt Pre→Post (ms)",
                ylabel="ΔW (%)", title="STDP curve")
                CairoMakie.hlines!(a4, [0.0]; color=:gray, linestyle=:dash)
                CairoMakie.vlines!(a4, [0.0]; color=:gray, linestyle=:dash)
                CairoMakie.scatter!(a4, sp.dt, sp.dw; color=:steelblue, markersize=10)
                CairoMakie.lines!(a4, sp.dt, sp.dw; color=:steelblue, linewidth=1.5)
                CairoMakie.save("plot_4_stdp.png", f4; px_per_unit=2)
                println("4 files exported (incl. STDP)")
            else
                println("3 files exported: plot_1_enzymes | plot_2_phase | plot_3_weight")
            end
        end

        # ── RUN callback ──────────────────────────────────────────────────────
        # Dispatches to the correct simulation based on active_tab[].
        on(session, btn_run.value) do _
            error_msg[] = ""
            tab = active_tab[]

            # ── Custom Stimulation ─────────────────────────────────────────────────
            if tab == "stdp_manual"
                v_ca,  err = get_valid_val(tb_man_ca,   2.0, 4.0,    "Ca²⁺");  err != "" && (error_msg[]=err; return)
                v_mg,  err = get_valid_val(tb_man_mg,   0.8, 3.0,    "Mg²⁺");  err != "" && (error_msg[]=err; return)
                v_tmp, err = get_valid_val(tb_man_temp, 20., 40.,    "Temp");   err != "" && (error_msg[]=err; return)
                v_pls, err = get_valid_val(tb_man_pls,  1.,  10000., "Pulses"); err != "" && (error_msg[]=err; return)
                v_frq, err = get_valid_val(tb_man_freq, 0.1, 100.,   "Freq");   err != "" && (error_msg[]=err; return)
                pre_ts,  err = parse_floats(tb_man_pre.value[],  "Pre times";  allow_negative=false)
                err != "" && (error_msg[]=err; return)
                post_ts, err = parse_floats(tb_man_post.value[], "Post times"; allow_negative=false)
                err != "" && (error_msg[]=err; return)
                isempty(pre_ts) && isempty(post_ts) && (error_msg[]="⚠ Pre and Post cannot both be empty."; return)
                v_nruns_man, err = get_valid_val(tb_n_runs, 1., 100., "N runs"); err != "" && (error_msg[]=err; return)
                println("▶ Custom — pre=$pre_ts ms  post=$post_ts ms  pulse=$(Int(v_pls))  freq=$v_frq Hz")

                current_task[] = errormonitor(@async begin
                    try
                        sim_fn = () -> begin
                            run_simulation!(p_bar, data, dataProtocol("TigaretMellor16"), 1;
                                plasticity_obs, use_protocol=false,
                                ca_ext=v_ca, mg=v_mg, temp=v_tmp,
                                pre_times=pre_ts, post_times=post_ts,
                                pulse=Int(v_pls), freq=v_frq, age=menu_man_age.value[])
                            return plasticity_obs[].weight[end]
                        end
                        run_batch!(Int(v_nruns_man), p_bar, plasticity_obs, final_wc, batch_results, all_run_curves, sim_fn)
                    catch e; e isa InterruptException || rethrow(e); end
                end)

            # ── HFS / LFS ──────────────────────────────────────────────────────
            elseif tab == "hfs"
                v_frq, err = get_valid_val(tb_hfs_freq,   0.1, 500.,    "Frequency"); err != "" && (error_msg[]=err; return)
                v_pls, err = get_valid_val(tb_hfs_pulses, 1.,  5000.,   "Pulses");    err != "" && (error_msg[]=err; return)
                v_ntr, err = get_valid_val(tb_hfs_ntr,    1.,  100.,    "N trains");  err != "" && (error_msg[]=err; return)
                v_itr, err = get_valid_val(tb_hfs_itr,    0.,  1000000., "Inter-train"); err != "" && (error_msg[]=err; return)
                v_ca,  err = get_valid_val(tb_hfs_ca,     2.0, 4.0,    "Ca2+");      err != "" && (error_msg[]=err; return)
                v_mg,  err = get_valid_val(tb_hfs_mg,     0.8, 3.0,    "Mg2+");      err != "" && (error_msg[]=err; return)
                v_tmp, err = get_valid_val(tb_hfs_temp,   20., 40.,    "Temp");      err != "" && (error_msg[]=err; return)
                v_nruns, err = get_valid_val(tb_n_runs, 1., 100., "N runs"); err != "" && (error_msg[]=err; return)


                mode_str    = v_frq >= 50 ? "HFS" : "LFS"
                n_trains    = Int(v_ntr)
                n_pulses    = Int(v_pls)
                inter_train = v_itr
                dt_pulse    = 1000.0 / v_frq   # ms between pulses within one train

                # Duration of one train: from first to last pulse
                train_duration = (n_pulses - 1) * dt_pulse

                println("$mode_str run: $n_trains trains x $n_pulses pulses @ $(v_frq) Hz, inter-train=$(inter_train) ms")

                # ── RAM safety check ──────────────────────────────────────────────────
                n_pre_hfs = n_trains * n_pulses
                t_end_hfs = n_trains * (train_duration + inter_train) + 1.5e5
                ok, ram_warn = check_memory_safety(n_pre_hfs, t_end_hfs)
                if !ok
                    error_msg[] = ram_warn
                    return
                end

                current_task[] = errormonitor(@async begin
                    try
                        start = 0.5e3

                        sim_fn = () -> begin
                            if n_trains == 1
                                run_simulation!(p_bar, data, dataProtocol("TigaretMellor16"), 1;
                                    plasticity_obs = plasticity_obs,
                                    use_protocol   = false,
                                    ca_ext = v_ca, mg = v_mg, temp = v_tmp,
                                    pre_times = [0.0], post_times = Float64[],
                                    pulse = n_pulses, freq = v_frq,
                                    age = menu_hfs_age.value[])
                            else
                                events_list = Tuple{Float64, Bool}[]
                                for tr in 0:(n_trains - 1)
                                    train_start = start + tr * (train_duration + inter_train)
                                    for p in 0:(n_pulses - 1)
                                        push!(events_list, (train_start + p * dt_pulse, true))
                                    end
                                end
                                sort!(events_list; by = x -> x[1])
                                events_times         = Float64[x[1] for x in events_list]
                                is_pre_or_post_event = Bool[x[2]   for x in events_list]
                                _ca_ext  = v_ca * 1000.0
                                _age     = (menu_hfs_age.value[] == "Adult" ? 60.0 : 21.0)
                                param_synapse = SynapseParams(
                                    t_end = events_times[end] + 1.5e5, soma_dist = 200.0,
                                    temp_rates = v_tmp, Ca_ext = _ca_ext, Mg = v_mg, age = _age,
                                    injbap = 2.0, I_clamp = 1000.0, sampling_rate = 0.01)
                                p_rel       = param_synapse.p_release
                                pre_synapse = PreSynapseParams(
                                    h = p_rel[1] + p_rel[2] / (1 + exp(p_rel[3] * (_ca_ext - p_rel[4]))))
                                xc0 = initial_conditions_continuous_temp(param_synapse)
                                xd0 = initial_conditions_discrete(param_synapse)
                                is_glu_release, _, _, _, _, bap_by_epsp_times =
                                    stp(param_synapse.t_end, pre_synapse, events_times,
                                        is_pre_or_post_event, _plot = false, algo = CHV(CVODE_BDF()))
                                valid_idx = events_times .< param_synapse.t_end
                                _evolve_synapse(
                                    xc0, xd0, param_synapse,
                                    events_times[valid_idx], is_pre_or_post_event[valid_idx],
                                    bap_by_epsp_times, is_glu_release,
                                    (CHV(CVODE_BDF(linear_solver = :GMRES)),
                                     CHV(CVODE_BDF(linear_solver = :GMRES)));
                                    abstol = 1e-6, reltol = 1e-5,
                                    save_positions = ((false, false), (false, true)),
                                    verbose = false, progress = p_bar,
                                    plasticity_obs = plasticity_obs)
                            end
                            return plasticity_obs[].weight[end]
                        end

                        run_batch!(Int(v_nruns), p_bar, plasticity_obs, final_wc, batch_results, all_run_curves, sim_fn)

                    catch e
                        e isa InterruptException || rethrow(e)
                    end
                end)


                # ── TBS ────────────────────────────────────────────────────────────
            elseif tab == "tbs"
                v_nb,  err = get_valid_val(tb_tbs_nb,  1.,  500.,   "N bursts");        err != "" && (error_msg[]=err; return)
                v_ppb, err = get_valid_val(tb_tbs_ppb, 1.,  20.,    "Pulses/burst");    err != "" && (error_msg[]=err; return)
                v_fi,  err = get_valid_val(tb_tbs_fi,  1.,  1000.,  "Intra frequency"); err != "" && (error_msg[]=err; return)
                v_ft,  err = get_valid_val(tb_tbs_ft,  0.1, 20.,    "Theta frequency"); err != "" && (error_msg[]=err; return)
                v_ntr, err = get_valid_val(tb_tbs_ntr, 1.,  100.,   "N trains");        err != "" && (error_msg[]=err; return)
                v_itr, err = get_valid_val(tb_tbs_itr, 0.,  1000000.,"Inter-train");     err != "" && (error_msg[]=err; return)
                v_ca,  err = get_valid_val(tb_tbs_ca,  2.0, 4.0,   "Ca²⁺");            err != "" && (error_msg[]=err; return)
                v_mg,  err = get_valid_val(tb_tbs_mg,  0.8, 3.0,   "Mg²⁺");            err != "" && (error_msg[]=err; return)
                v_tmp, err = get_valid_val(tb_tbs_temp, 20., 40.,  "Temp");             err != "" && (error_msg[]=err; return)

                # Durées élémentaires
                dt_intra    = 1000.0 / v_fi          # ms entre pulses dans un burst
                dt_theta    = 1000.0 / v_ft          # ms entre bursts dans un train
                n_bursts    = Int(v_nb)
                n_pulses    = Int(v_ppb)
                n_trains    = Int(v_ntr)
                inter_train = v_itr                  # ms entre les trains

                # Construction des offsets d'un burst (relatifs à son début)
                burst_offsets = [p * dt_intra for p in 0:(n_pulses - 1)]
                # Durée d'un train entier
                train_duration = (n_bursts - 1) * dt_theta + burst_offsets[end]

                println(" TBS — $n_trains trains × $n_bursts bursts × $n_pulses pulses")
                println("  Intra: $(v_fi) Hz | Theta: $(v_ft) Hz | Inter-train: $(inter_train) ms")
                println("  Burst offsets: $burst_offsets ms")
                println("  Train duration: $(round(train_duration; digits=1)) ms")

                # ── RAM safety check ──────────────────────────────────────────────────
                n_pre_tbs = n_trains * n_bursts * n_pulses
                t_end_tbs = n_trains * (train_duration + inter_train) + 1.5e5
                ok, ram_warn = check_memory_safety(n_pre_tbs, t_end_tbs)
                if !ok
                    error_msg[] = ram_warn
                    return
                end

                current_task[] = errormonitor(@async begin
                    try
                        start = 0.5e3   # 500 ms offset initial

                        # ── Construction manuelle de TOUS les événements ──────────────────
                        # Pour chaque train, pour chaque burst, pour chaque pulse :
                        # on calcule le temps absolu et on l'ajoute à la liste.
                        events_list = Tuple{Float64, Bool}[]

                        for tr in 0:(n_trains - 1)
                            # Début de ce train (ms)
                            train_start = start + tr * (train_duration + inter_train)

                            for b in 0:(n_bursts - 1)
                                # Début de ce burst dans le train
                                burst_start = train_start + b * dt_theta

                                for offset in burst_offsets
                                    # Chaque pulse est un événement pré-synaptique
                                    push!(events_list, (burst_start + offset, true))
                                end
                            end
                        end

                        # Tri chronologique obligatoire pour le solveur PDMP
                        sort!(events_list; by = x -> x[1])

                        events_times         = Float64[x[1] for x in events_list]
                        is_pre_or_post_event = Bool[x[2]   for x in events_list]

                        println("  Total events: $(length(events_times))")
                        println("  First event: $(events_times[1]) ms")
                        println("  Last event:  $(events_times[end]) ms")

                        # ── Paramètres synaptiques ────────────────────────────────────────
                        _ca_ext = v_ca * 1000.0
                        _age    = (menu_tbs_age.value[] == "Adult" ? 60.0 : 21.0)
                        afterend = 1.5e5   # 150 s après le dernier spike

                        param_synapse = SynapseParams(
                        t_end      = events_times[end] + afterend,
                        soma_dist  = 200.0,
                        temp_rates = v_tmp,
                        Ca_ext     = _ca_ext,
                        Mg         = v_mg,
                        age        = _age,
                        injbap     = 2.0,
                        I_clamp    = 1000.0,
                        sampling_rate = 0.01)

                        p           = param_synapse.p_release
                        pre_synapse = PreSynapseParams(
                        h = p[1] + p[2] / (1 + exp(p[3] * (_ca_ext - p[4]))))

                        xc0 = initial_conditions_continuous_temp(param_synapse)
                        xd0 = initial_conditions_discrete(param_synapse)

                        # Short-term plasticity pré-synaptique
                        is_glu_release, _, _, _, _, bap_by_epsp_times =
                        stp(param_synapse.t_end, pre_synapse,
                        events_times, is_pre_or_post_event,
                        _plot = false, algo = CHV(CVODE_BDF()))

                        valid_idx = events_times .< param_synapse.t_end
                        ev_clean  = events_times[valid_idx]
                        ipp_clean = is_pre_or_post_event[valid_idx]

                        # ── Simulation ────────────────────────────────────────────────────
                        v_nruns_tbs = tryparse(Float64, tb_n_runs.value[])
                        v_nruns_tbs = isnothing(v_nruns_tbs) ? 1 : Int(clamp(v_nruns_tbs, 1., 100.))

                        sim_fn = () -> begin
                            _evolve_synapse(
                                xc0, xd0, param_synapse, ev_clean, ipp_clean,
                                bap_by_epsp_times, is_glu_release,
                                (CHV(CVODE_BDF(linear_solver = :GMRES)),
                                 CHV(CVODE_BDF(linear_solver = :GMRES)));
                                abstol = 1e-6, reltol = 1e-5,
                                save_positions = ((false, false), (false, true)),
                                verbose = false, progress = p_bar,
                                plasticity_obs = plasticity_obs)
                            return plasticity_obs[].weight[end]
                        end

                        run_batch!(v_nruns_tbs, p_bar, plasticity_obs, final_wc, batch_results, all_run_curves, sim_fn)
                        @warn "TBS finished"
                    catch e
                        e isa InterruptException || rethrow(e)
                    end
                end)

            # ── Protocol ───────────────────────────────────────────────────────
            elseif tab == "protocol"
                k  = parse(Int, string(menu_exp.value[]))
                dp = dataProtocol(menu_prot.value[])
                println("▶ Protocol: $(menu_prot.value[])  experiment #$k")

                # ── RAM safety check ──────────────────────────────────────────────────
                # Read simulation parameters directly from the protocol DataFrame row
                _pls_ram   = dp[!, :pulse][k]
                _freq_ram  = dp[!, :freq][k]
                _rept_ram  = dp[!, :repeat_times][k]
                n_pre_prot = (_rept_ram + 1) * _pls_ram          # total pre events
                t_end_prot = 0.5e3 + 0.1 * (_rept_ram + 1) * _pls_ram * 1000.0 / _freq_ram + 1.5e5
                ok, ram_warn = check_memory_safety(n_pre_prot, t_end_prot)
                if !ok
                    error_msg[] = ram_warn
                    return
                end

                v_nruns_prot, err = get_valid_val(tb_n_runs, 1., 100., "N runs"); err != "" && (error_msg[]=err; return)
                current_task[] = errormonitor(@async begin
                    try
                        sim_fn = () -> begin
                            run_simulation!(p_bar, data, dp, k;
                                plasticity_obs, use_protocol=true)
                            return plasticity_obs[].weight[end]
                        end
                        run_batch!(Int(v_nruns_prot), p_bar, plasticity_obs, final_wc, batch_results, all_run_curves, sim_fn)
                    catch e; e isa InterruptException || rethrow(e); end
                end)
            end

        end   # end on(btn_run)

                # ── DOM ───────────────────────────────────────────────────────────────

                # Tab bar styling constants
                TAB_BASE   = "padding:9px 18px;border:none;border-radius:8px 8px 0 0;" *
                "font-size:13px;font-weight:600;cursor:pointer;" *
                "font-family:'Segoe UI',sans-serif;transition:all 0.15s;"
                TAB_ACTIVE = TAB_BASE * "background:#2c3e50;color:white;"
                TAB_IDLE   = TAB_BASE * "background:#ecf0f1;color:#7f8c8d;"

                # The tab bar: 5 buttons side by side.
                # We cannot change button styles reactively from Julia in Bonito,
                # so all buttons use the same neutral style — the active tab is
                # indicated by the content area below, not by button colour.
                tab_bar = Bonito.DOM.div(
                    btn_tab_hfs, btn_tab_tbs, btn_tab_stdp_manual, btn_tab_protocol;
                    style = "display:flex;gap:4px;padding:0 0 0 2px;" *
                            "border-bottom:3px solid #2c3e50;margin-bottom:0;")

                # Tab content area: lift(active_tab) re-evaluates whenever the tab changes.
                # Each branch calls the corresponding tab_xxx() builder function.
                # This is the key mechanism: the DOM node returned changes, Bonito
                # re-renders only this subtree. All widgets are already registered
                # (they were instantiated above), so Bonito never crashes looking for
                # a missing data-jscall-id.
                tab_content = lift(active_tab) do tab
                    if     tab == "hfs"        ; tab_hfs()
                    elseif tab == "tbs"        ; tab_tbs()
                    elseif tab == "stdp_manual"; tab_stdp_manual()
                    elseif tab == "protocol"   ; tab_protocol()
                    else                       ; Bonito.DOM.div()
                    end
                end


                return Bonito.DOM.div(

                # ── Page header ───────────────────────────────────────────────────
                Bonito.DOM.h1("SynapseElife Studio";
                style = "color:#2c3e50;font-family:'Segoe UI',sans-serif;" *
                "text-align:center;margin-bottom:24px;letter-spacing:1px;"),

               # ── Global toolbar + final weight change badge ────────────────────────────
                Bonito.DOM.div(
                    # Left side: action buttons
                    Bonito.DOM.div(
                        btn_run, btn_stop, btn_era,
                        Bonito.DOM.span(" | "; style="color:#ccc;margin:0 8px;"),
                        btn_csv, btn_plt,
                        # N runs field inline in the toolbar
                        Bonito.DOM.span("N runs:"; style="font-size:12px;color:#7f8c8d;margin-left:12px;"),
                        tb_n_runs;
                        style = "display:flex;gap:8px;align-items:center;"),
                    # Right side: final weight change badge — updated at end of each simulation
                    Bonito.DOM.div(final_wc;
                        style = "font-size:15px;font-weight:bold;color:#2c3e50;" *
                        "background:#eaf4fb;padding:6px 16px;border-radius:8px;" *
                        "border:1px solid #aed6f1;white-space:nowrap;");
                    # Outer toolbar div style — flex row with space-between
                        style = "background:white;padding:14px 20px;border-radius:10px;" *
                        "box-shadow:0 3px 10px rgba(0,0,0,0.09);margin-bottom:14px;" *
                        "display:flex;justify-content:space-between;align-items:center;"),
                
                # Status messages
                Bonito.DOM.div(confirm_stop_msg;
                style = "color:#d97706;font-weight:600;font-size:13px;" *
                "text-align:center;min-height:1.2em;margin-bottom:4px;"),
                Bonito.DOM.div(error_msg;
                style = "color:#e74c3c;font-weight:600;font-size:13px;" *
                "text-align:center;min-height:1.2em;margin-bottom:10px;"),
                
                # ── Tab system ────────────────────────────────────────────────────
                Bonito.DOM.div(
                tab_bar,
                # Content panel — updates reactively when active_tab changes
                Bonito.DOM.div(tab_content;
                style = "background:white;border-radius:0 10px 10px 10px;" *
                "box-shadow:0 3px 10px rgba(0,0,0,0.09);");
                style = "margin-bottom:16px;"),
                
                # ── Main simulation figure ─────────────────────────────────────────
                Bonito.DOM.div(fig;
                style = "background:white;padding:16px;border-radius:10px;" *
                "box-shadow:0 3px 10px rgba(0,0,0,0.09);" *
                "display:flex;justify-content:center;");
                
                style = "background-color:#f0f2f5;padding:32px;" *
                "min-height:100vh;box-sizing:border-box;" *
                "font-family:'Segoe UI',sans-serif;max-width:1300px;margin:0 auto;"
                )
            end   # end Bonito.App() do session
            return app
        end
        
        
        # =============================================================================
        #  LAUNCH
        # =============================================================================
        app    = build_app()
        server = Bonito.Server(app, "127.0.0.1", 8081)
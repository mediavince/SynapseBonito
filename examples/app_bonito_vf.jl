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
#    • Manual   : fully custom pre/post spike timestamps
#
#  EXPORT
#  ──────
#  WGLMakie renders in the browser via WebGL — it cannot write PNG files.
#  CairoMakie (server-side) is used exclusively for file export.
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
Statistics

import CairoMakie

includet("utils_bonito.jl")


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

# Render a protocol DataFrame as a scrollable HTML table.
# The row at selected_idx is highlighted in blue.
function format_table(dp, selected_idx = 1)
    df = copy(dp)
    DataFrames.insertcols!(df, 1, "ID" => 1:DataFrames.nrow(df))

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
        stdp_points      = Observable((dt = Float64[], dw = Float64[]))

        # active_tab: the single source of truth for which tab is shown.
        # Possible values: "stdp" | "hfs" | "tbs" | "manual" | "protocol"
        # Default is "stdp".
        active_tab = Observable("stdp")

        # ── GLOBAL TOOLBAR WIDGETS ────────────────────────────────────────────
        btn_run  = Bonito.Button("▶  Run")
        btn_stop = Bonito.Button("■  Stop")
        btn_era  = Bonito.Button("✕  Erase")
        btn_csv  = Bonito.Button("Export CSV")
        btn_plt  = Bonito.Button("Export Plots")

        # ── TAB BUTTONS ───────────────────────────────────────────────────────
        # Each button sets active_tab[] to the corresponding tab name.
        # They are plain Bonito.Button widgets — no JS onclick needed.
        btn_tab_stdp     = Bonito.Button("STDP")
        btn_tab_hfs      = Bonito.Button("HFS / LFS")
        btn_tab_tbs      = Bonito.Button("TBS")
        btn_tab_manual   = Bonito.Button("Manual")
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
        tb_hfs_ca     = Bonito.TextField("2.5")
        tb_hfs_mg     = Bonito.TextField("1.0")
        tb_hfs_temp   = Bonito.TextField("35.0")
        menu_hfs_age  = Bonito.Dropdown(["Young", "Adult"]; index = 2)

        # ── STDP TAB WIDGETS ──────────────────────────────────────────────────
        # Δt convention: Δt = t_post − t_pre
        #   Δt > 0 → Pre before Post → causal → LTP expected
        #   Δt < 0 → Post before Pre → anti-causal → LTD expected
        # One independent simulation is run per Δt value; the final weight
        # of each run is one data point on the STDP curve.
        tb_stdp_dt     = Bonito.TextField("-50 -20 -10 -5 5 10 20 50")
        tb_stdp_rep    = Bonito.TextField("60")      # repetitions per Δt
        tb_stdp_interv = Bonito.TextField("1000.0")  # inter-burst interval (ms)
        tb_stdp_ca     = Bonito.TextField("2.5")
        tb_stdp_mg     = Bonito.TextField("1.0")
        tb_stdp_temp   = Bonito.TextField("35.0")
        menu_stdp_age  = Bonito.Dropdown(["Young", "Adult"]; index = 2)
        btn_stdp_clear = Bonito.Button("Clear curve")

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

        lines!(ax_wt, lift(x -> x.t./1000, plasticity_obs),
        lift(x -> x.weight,  plasticity_obs))

        poly!(ax_phase, Point2f[(P[1],P[2]) for P in SynapseParams().LTP_region.vertices];
        color = ("red",  0.1))
        poly!(ax_phase, Point2f[(P[1],P[2]) for P in SynapseParams().LTD_region.vertices];
        color = ("blue", 0.1))

        # STDP figure — embedded in the STDP tab content
        fig_stdp = Figure(; size = (480, 320), fontsize = fontsize * 0.8)
        ax_stdp  = Axis(fig_stdp[1, 1];
        xlabel = "Δt  Pre→Post  (ms)", ylabel = "ΔW (%)", title = "STDP curve",
        ylabelsize = fontsize*0.55, xlabelsize = fontsize*0.55,
        yticklabelsize = fontsize*0.5, xticklabelsize = fontsize*0.5)
        hlines!(ax_stdp, [0.0]; color = :gray, linestyle = :dash, linewidth = 1)
        vlines!(ax_stdp, [0.0]; color = :gray, linestyle = :dash, linewidth = 1)
        scatter!(ax_stdp, lift(x -> x.dt, stdp_points),
        lift(x -> x.dw, stdp_points);
        color = :steelblue, markersize = 10)
        lines!(ax_stdp,   lift(x -> x.dt, stdp_points),
        lift(x -> x.dw, stdp_points);
        color = :steelblue, linewidth = 1.5)

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
        HINT   = "font-size:11px;color:#95a5a6;margin-bottom:12px;line-height:1.5;"

        # Helper: 4 shared physical-parameter fields used by all stimulation tabs
        function phys_fields(tb_c, tb_m, tb_t, menu_a)
            Bonito.DOM.div(
            Bonito.DOM.div(Bonito.DOM.span("Ca²⁺ [2.5–4] (mM):", style=LBL), tb_c,   style=COL),
            Bonito.DOM.div(Bonito.DOM.span("Mg²⁺ [0.8–3] (mM):", style=LBL), tb_m,   style=COL),
            Bonito.DOM.div(Bonito.DOM.span("Temp [20–40] (°C):", style=LBL),  tb_t,   style=COL),
            Bonito.DOM.div(Bonito.DOM.span("Age:",               style=LBL),  menu_a, style=COL);
            style = GRID)
        end

        # ── TAB CONTENT BUILDERS ──────────────────────────────────────────────
        # Each function returns a DOM node for that tab's parameter panel.
        # They are called inside lift(active_tab) so only the active one
        # is rendered in the visible part of the page.

        function tab_stdp()
            Bonito.DOM.div(
            Bonito.DOM.h3("Spike Timing-Dependent Plasticity"; style=TITLE),
            Bonito.DOM.p(
            "Enter Δt values (ms) separated by spaces. " *
            "Δt > 0 → Pre fires before Post (causal, LTP expected). " *
            "Δt < 0 → Post fires before Pre (anti-causal, LTD expected). " *
            "One independent simulation runs per Δt; the final weight " *
            "of each run becomes one point on the STDP curve below.";
            style = HINT),
            Bonito.DOM.div(
            Bonito.DOM.div(Bonito.DOM.span("Δt list (ms):",   style=LBL*"font-weight:bold;"), tb_stdp_dt,    style=COL),
            Bonito.DOM.div(Bonito.DOM.span("Repetitions:",    style=LBL), tb_stdp_rep,    style=COL),
            Bonito.DOM.div(Bonito.DOM.span("Interval (ms):",  style=LBL), tb_stdp_interv, style=COL),
            phys_fields(tb_stdp_ca, tb_stdp_mg, tb_stdp_temp, menu_stdp_age);
            style=GRID),
            Bonito.DOM.div(btn_stdp_clear; style="margin-top:12px;"),
            Bonito.DOM.div(fig_stdp; style="margin-top:18px;display:flex;justify-content:center;");
            style = CARD)
        end

        function tab_hfs()
            Bonito.DOM.div(
            Bonito.DOM.h3("High / Low Frequency Stimulation"; style=TITLE),
            Bonito.DOM.p(
            "Simulates a pure pre-synaptic spike train (no post-synaptic spikes). " *
            "Frequency ≥ 50 Hz is classified as HFS; < 50 Hz as LFS.";
            style = HINT),
            Bonito.DOM.div(
            Bonito.DOM.div(Bonito.DOM.span("Frequency (Hz):", style=LBL), tb_hfs_freq,   style=COL),
            Bonito.DOM.div(Bonito.DOM.span("Nb pulses:",      style=LBL), tb_hfs_pulses, style=COL),
            phys_fields(tb_hfs_ca, tb_hfs_mg, tb_hfs_temp, menu_hfs_age);
            style=GRID);
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
            Bonito.DOM.div(
            Bonito.DOM.div(Bonito.DOM.span("N trains:",            style=LBL), tb_tbs_ntr, style=COL),
            Bonito.DOM.div(Bonito.DOM.span("Inter-train (ms):",    style=LBL), tb_tbs_itr, style=COL),
            Bonito.DOM.div(Bonito.DOM.span("N bursts / train:",    style=LBL), tb_tbs_nb,  style=COL),
            Bonito.DOM.div(Bonito.DOM.span("Pulses / burst:",      style=LBL), tb_tbs_ppb, style=COL),
            Bonito.DOM.div(Bonito.DOM.span("Intra freq (Hz):",     style=LBL), tb_tbs_fi,  style=COL),
            Bonito.DOM.div(Bonito.DOM.span("Theta freq (Hz):",     style=LBL), tb_tbs_ft,  style=COL),
            phys_fields(tb_tbs_ca, tb_tbs_mg, tb_tbs_temp, menu_tbs_age);
            style=GRID);
            style = CARD * "border-left:4px solid #b45309;")
        end

        function tab_manual()
            Bonito.DOM.div(
            Bonito.DOM.h3("Manual / Custom Pattern"; style=TITLE),
            Bonito.DOM.p(
            "Enter absolute spike timestamps (ms) within one burst cycle, " *
            "separated by spaces. The pattern is repeated `Pulses` times " *
            "at the given frequency. " *
            "Example: Pre = \"0 10 20\", Post = \"30\" → 3 pre spikes + 1 post spike.";
            style = HINT),
            Bonito.DOM.div(
            Bonito.DOM.div(Bonito.DOM.span("Pre times (ms):",     style=LBL*"font-weight:bold;"), tb_man_pre,  style=COL),
            Bonito.DOM.div(Bonito.DOM.span("Post times (ms):",    style=LBL*"font-weight:bold;"), tb_man_post, style=COL),
            Bonito.DOM.div(Bonito.DOM.span("Repetitions:",        style=LBL), tb_man_pls,  style=COL),
            Bonito.DOM.div(Bonito.DOM.span("Frequency (Hz):",     style=LBL), tb_man_freq, style=COL),
            phys_fields(tb_man_ca, tb_man_mg, tb_man_temp, menu_man_age);
            style=GRID);
            style = CARD * "border-left:4px solid #c0392b;")
        end

        function tab_protocol()
            Bonito.DOM.div(
            Bonito.DOM.h3("Protocol Mode"; style=TITLE),
            Bonito.DOM.p(
            "Select a published protocol from the dropdown. " *
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
        on(session, btn_tab_stdp.value)     do _; active_tab[] = "stdp";     end
        on(session, btn_tab_hfs.value)      do _; active_tab[] = "hfs";      end
        on(session, btn_tab_tbs.value)      do _; active_tab[] = "tbs";      end
        on(session, btn_tab_manual.value)   do _; active_tab[] = "manual";   end
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

        on(session, btn_stdp_clear.value) do _
            stdp_points[] = (dt = Float64[], dw = Float64[])
        end

        # Export CSV
        on(session, btn_csv.value) do _
            obs = plasticity_obs[]
            CSV.write("resultats_synapse.csv", DataFrame(
            Time_s               = obs.t ./ 1000,
            CaMKII_uM            = obs.camk2,
            CaN_uM               = obs.can,
            WeightChange_percent = obs.weight))
            if !isempty(stdp_points[].dt)
                CSV.write("resultats_stdp.csv", DataFrame(
                Delta_t_ms       = stdp_points[].dt,
                Weight_change_pc = stdp_points[].dw))
                println("CSV → resultats_synapse.csv + resultats_stdp.csv")
            else
                println("CSV → resultats_synapse.csv")
            end
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

            # ── STDP ───────────────────────────────────────────────────────────
            if tab == "stdp"
                dt_list, err = parse_floats(tb_stdp_dt.value[], "Δt"; allow_negative=true)
                err != "" && (error_msg[]=err; return)
                isempty(dt_list) && (error_msg[]="⚠ Enter at least one Δt value."; return)

                v_rep, err = get_valid_val(tb_stdp_rep,    1.,   1000.,   "Repetitions"); err != "" && (error_msg[]=err; return)
                v_int, err = get_valid_val(tb_stdp_interv, 10.,  100000., "Interval");    err != "" && (error_msg[]=err; return)
                v_ca,  err = get_valid_val(tb_stdp_ca,  2.5, 4.0, "Ca²⁺"); err != "" && (error_msg[]=err; return)
                v_mg,  err = get_valid_val(tb_stdp_mg,  0.8, 3.0, "Mg²⁺"); err != "" && (error_msg[]=err; return)
                v_tmp, err = get_valid_val(tb_stdp_temp, 20., 40., "Temp"); err != "" && (error_msg[]=err; return)

                sort!(dt_list)
                println("▶ STDP — $(length(dt_list)) Δt values: $dt_list ms")

                current_task[] = errormonitor(@async begin
                    try
                        for dt in dt_list
                            # Reset plots before each individual Δt simulation
                            plasticity_obs[] = (t=zeros(1), camk2=zeros(1),
                            can=zeros(1), weight=zeros(1))
                            p_bar[] = 0.0

                            # Δt ≥ 0 → pre at 0 ms, post at Δt ms (causal)
                            # Δt < 0 → post at 0 ms, pre at |Δt| ms (anti-causal)
                            pre_ts, post_ts = dt >= 0 ? ([0.0], [dt]) : ([abs(dt)], [0.0])

                            w = run_one_simulation!(p_bar, plasticity_obs;
                            pre_ts=pre_ts, post_ts=post_ts,
                            repeat_n=Int(v_rep), repeat_interval=v_int,
                            ca_ext=v_ca, mg=v_mg, temp=v_tmp,
                            age=menu_stdp_age.value[])

                            # Append point and keep sorted by Δt
                            old = stdp_points[]
                            ndt = vcat(old.dt, dt); ndw = vcat(old.dw, w)
                            pm  = sortperm(ndt)
                            stdp_points[] = (dt=ndt[pm], dw=ndw[pm])
                            println("  Δt=$(dt) ms → ΔW=$(round(w; digits=2))%")
                        end
                        println("STDP curve complete!")
                    catch e; e isa InterruptException || rethrow(e); end
                end)

                # ── HFS / LFS ──────────────────────────────────────────────────────
            elseif tab == "hfs"
                v_frq, err = get_valid_val(tb_hfs_freq,   0.1, 500.,  "Frequency"); err != "" && (error_msg[]=err; return)
                v_pls, err = get_valid_val(tb_hfs_pulses, 1.,  5000., "Pulses");    err != "" && (error_msg[]=err; return)
                v_ca,  err = get_valid_val(tb_hfs_ca,  2.5, 4.0, "Ca²⁺"); err != "" && (error_msg[]=err; return)
                v_mg,  err = get_valid_val(tb_hfs_mg,  0.8, 3.0, "Mg²⁺"); err != "" && (error_msg[]=err; return)
                v_tmp, err = get_valid_val(tb_hfs_temp, 20., 40., "Temp"); err != "" && (error_msg[]=err; return)

                mode_str = v_frq >= 50 ? "HFS" : "LFS"
                println("▶ $mode_str — freq=$(v_frq) Hz  pulses=$(Int(v_pls))")
                current_task[] = errormonitor(@async begin
                    try
                        run_simulation!(p_bar, data, dataProtocol("TigaretMellor16"), 1;
                        plasticity_obs, use_protocol=false,
                        ca_ext=v_ca, mg=v_mg, temp=v_tmp,
                        pre_times=[0.0], post_times=Float64[],
                        pulse=Int(v_pls), freq=v_frq, age=menu_hfs_age.value[])
                    catch e; e isa InterruptException || rethrow(e); end
                end)

                # ── TBS ────────────────────────────────────────────────────────────
            elseif tab == "tbs"
                v_nb,  err = get_valid_val(tb_tbs_nb,  1.,  500.,   "N bursts");        err != "" && (error_msg[]=err; return)
                v_ppb, err = get_valid_val(tb_tbs_ppb, 1.,  20.,    "Pulses/burst");    err != "" && (error_msg[]=err; return)
                v_fi,  err = get_valid_val(tb_tbs_fi,  1.,  1000.,  "Intra frequency"); err != "" && (error_msg[]=err; return)
                v_ft,  err = get_valid_val(tb_tbs_ft,  0.1, 20.,    "Theta frequency"); err != "" && (error_msg[]=err; return)
                v_ntr, err = get_valid_val(tb_tbs_ntr, 1.,  100.,   "N trains");        err != "" && (error_msg[]=err; return)
                v_itr, err = get_valid_val(tb_tbs_itr, 0.,  300000.,"Inter-train");     err != "" && (error_msg[]=err; return)
                v_ca,  err = get_valid_val(tb_tbs_ca,  2.5, 4.0,   "Ca²⁺");            err != "" && (error_msg[]=err; return)
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
                        @time _evolve_synapse(
                        xc0, xd0, param_synapse, ev_clean, ipp_clean,
                        bap_by_epsp_times, is_glu_release,
                        (CHV(CVODE_BDF(linear_solver = :GMRES)),
                        CHV(CVODE_BDF(linear_solver = :GMRES)));
                        abstol = 1e-6, reltol = 1e-5,
                        save_positions = ((false, false), (false, true)),
                        verbose = false, progress = p_bar, plasticity_obs = plasticity_obs)

                        @warn "TBS finished"
                    catch e
                        e isa InterruptException || rethrow(e)
                    end
                end)

            # ── Manual ─────────────────────────────────────────────────────────
            elseif tab == "manual"
                v_ca,  err = get_valid_val(tb_man_ca,   2.5, 4.0,    "Ca²⁺");  err != "" && (error_msg[]=err; return)
                v_mg,  err = get_valid_val(tb_man_mg,   0.8, 3.0,    "Mg²⁺");  err != "" && (error_msg[]=err; return)
                v_tmp, err = get_valid_val(tb_man_temp, 20., 40.,    "Temp");   err != "" && (error_msg[]=err; return)
                v_pls, err = get_valid_val(tb_man_pls,  1.,  10000., "Pulses"); err != "" && (error_msg[]=err; return)
                v_frq, err = get_valid_val(tb_man_freq, 0.1, 100.,   "Freq");   err != "" && (error_msg[]=err; return)

                pre_ts,  err = parse_floats(tb_man_pre.value[],  "Pre times";  allow_negative=false)
                err != "" && (error_msg[]=err; return)
                post_ts, err = parse_floats(tb_man_post.value[], "Post times"; allow_negative=false)
                err != "" && (error_msg[]=err; return)
                isempty(pre_ts) && isempty(post_ts) &&
                    (error_msg[]="⚠ Pre and Post cannot both be empty."; return)

                println("▶ Manual — pre=$pre_ts ms  post=$post_ts ms  pulse=$(Int(v_pls))  freq=$v_frq Hz")
                current_task[] = errormonitor(@async begin
                    try
                        run_simulation!(p_bar, data, dataProtocol("TigaretMellor16"), 1;
                            plasticity_obs, use_protocol=false,
                            ca_ext=v_ca, mg=v_mg, temp=v_tmp,
                            pre_times=pre_ts, post_times=post_ts,
                            pulse=Int(v_pls), freq=v_frq, age=menu_man_age.value[])
                    catch e; e isa InterruptException || rethrow(e); end
                end)

            # ── Protocol ───────────────────────────────────────────────────────
            elseif tab == "protocol"
                k  = parse(Int, string(menu_exp.value[]))
                dp = dataProtocol(menu_prot.value[])
                println("▶ Protocol: $(menu_prot.value[])  experiment #$k")
                current_task[] = errormonitor(@async begin
                    try
                        run_simulation!(p_bar, data, dp, k;
                            plasticity_obs, use_protocol=true)
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
                btn_tab_stdp, btn_tab_hfs, btn_tab_tbs, btn_tab_manual, btn_tab_protocol;
                style = "display:flex;gap:4px;padding:0 0 0 2px;" *
                "border-bottom:3px solid #2c3e50;margin-bottom:0;")

                # Tab content area: lift(active_tab) re-evaluates whenever the tab changes.
                # Each branch calls the corresponding tab_xxx() builder function.
                # This is the key mechanism: the DOM node returned changes, Bonito
                # re-renders only this subtree. All widgets are already registered
                # (they were instantiated above), so Bonito never crashes looking for
                # a missing data-jscall-id.
                tab_content = lift(active_tab) do tab
                    if     tab == "stdp"     ; tab_stdp()
                    elseif tab == "hfs"      ; tab_hfs()
                    elseif tab == "tbs"      ; tab_tbs()
                    elseif tab == "manual"   ; tab_manual()
                    elseif tab == "protocol" ; tab_protocol()
                    else                     ; Bonito.DOM.div()
                    end
                end

                return Bonito.DOM.div(

                # ── Page header ───────────────────────────────────────────────────
                Bonito.DOM.h1("SynapseElife Studio";
                style = "color:#2c3e50;font-family:'Segoe UI',sans-serif;" *
                "text-align:center;margin-bottom:24px;letter-spacing:1px;"),

                # ── Global toolbar ────────────────────────────────────────────────
                Bonito.DOM.div(
                Bonito.DOM.div(
                btn_run, btn_stop, btn_era,
                Bonito.DOM.span(" | "; style="color:#ccc;margin:0 8px;"),
                btn_csv, btn_plt;
                style = "display:flex;gap:8px;align-items:center;"),
                style = "background:white;padding:14px 20px;border-radius:10px;" *
                "box-shadow:0 3px 10px rgba(0,0,0,0.09);margin-bottom:14px;"),
                
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
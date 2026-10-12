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
CSV

import CairoMakie

includet("utils_bonito.jl")

# ──────────────────────────────────────────────────────────────────
#  API Bonito (version stable 2024-2025)
#  • Dropdown(options)          → liste déroulante, .value :: Observable
#  • Button("label")            → bouton,           .value :: Observable{Bool}
#  • on(session, obs) do ...    → callback lié à la session (évite les fuites)
#  • DOM.div(...)               → nœud HTML
#  • display(App() do ... end)  → affiche dans VSCode / navigateur
# ──────────────────────────────────────────────────────────────────

function run_simulation!(progress, data, data_protocol, k;
                            plasticity_obs,
                            use_protocol=true,
                            ca_ext=2.5,
                            mg=1.0,
                            temp=35.0,
                            pre_times=[0.0],
                            post_times=[10.0],
                            pulse=300,
                            freq=5.0, age="Adult")

    _debug = true

    # ── DATA ROUTING & PATTERN GENERATION ──
    start = 0.5e3

    if use_protocol
        @info "Protocol mode"
        pls           = data_protocol[!,:pulse][k]
        _freq         = data_protocol[!,:freq][k]
        _temp         = data_protocol[!,:temp][k]
        _ca_ext       = data_protocol[!,:exca][k]
        _mg           = data_protocol[!,:exmg][k]
        _age          = data_protocol[!,:age][k]
        _repeat_times = data_protocol[!,:repeat_times][k]
        _injbap       = data_protocol[!,:inj_time][k]
        _injection    = data_protocol[!,:injection][k]
        _AP_by_EPSP   = data_protocol[!,:AP_by_EPSP][k]

        # Protocol Mode
        events_times, is_pre_or_post_event = firingPattern(
            start_time   = start,
            n_pos        = data_protocol[!,:n_pos][k],
            delay_pos    = data_protocol[!,:delay_pos][k],
            n_pre        = data_protocol[!,:n_pre][k],
            delay_pre    = data_protocol[!,:delay_pre][k],
            delay        = data_protocol[!,:delay][k],
            pulse        = pls,
            freq         = _freq,
            causal       = data_protocol[!,:causal][k],
            repeat_times = _repeat_times,
            repeat_after = data_protocol[!,:repeat_after][k]
            )

    else
        @info "Manual mode"
        pls           = pulse
        _freq         = freq
        _temp         = temp
        _ca_ext       = ca_ext * 1000.0
        _mg           = mg
        _age          = (age == "Adult" ? 60.0 : 21.0)
        _repeat_times = 0
        _injbap       = 2.0
        _injection    = 1000.0
        _AP_by_EPSP   = "yes"

        # ---IDEA ---
        events_list = Tuple{Float64, Bool}[]

        for p in 0:(pls - 1)
            base_t = start + p * (1000.0 / freq)

            # Add pre-synaptic events
            for t_pre in pre_times
                push!(events_list, (base_t + t_pre, true))
            end

            # Add post-synaptic events
            for t_post in post_times
                push!(events_list, (base_t + t_post, false))
            end
        end

        # Sort chronologically
        sort!(events_list, by = x -> x[1])

        # Split into 2 lists
        events_times = Float64[x[1] for x in events_list]
        is_pre_or_post_event = Bool[x[2] for x in events_list]
    end

    afterend      = 1.5e5
    param_synapse = SynapseParams(
    t_end         = start + 0.1 * (_repeat_times + 1) * pls * 1000 / _freq + afterend,
    soma_dist     = 200.0,
    temp_rates    = _temp,
    Ca_ext        = _ca_ext,
    Mg            = _mg,
    age           = _age,
    injbap        = _injbap,
    I_clamp       = _injection,
    sampling_rate = 0.01)

    p           = param_synapse.p_release

    pre_synapse = PreSynapseParams(
    h = p[1] + p[2] / (1 + exp(p[3] * (_ca_ext - p[4]))))

    xc0 = initial_conditions_continuous_temp(param_synapse)
    xd0 = initial_conditions_discrete(param_synapse)

    is_glu_release, Docked, Reserve, t_stp, glu_release_times, bap_by_epsp_times =
    stp(param_synapse.t_end, pre_synapse, events_times, is_pre_or_post_event,
    _plot = false, algo = CHV(CVODE_BDF()))

    @show "number of releases $(sum(is_glu_release))"

    # SYNAPSE EVOLUTION
    valid_idx = events_times .< param_synapse.t_end
    events_sorted_times_clean = events_times[valid_idx]
    is_pre_or_post_event_clean = is_pre_or_post_event[valid_idx]

    if _debug
        result = @time _evolve_synapse(
        xc0, xd0,
        param_synapse,
        events_sorted_times_clean,
        is_pre_or_post_event_clean,
        _AP_by_EPSP == "yes" ? bap_by_epsp_times : Float64[],
        is_glu_release,
        (CHV(CVODE_BDF(linear_solver = :GMRES)), CHV(CVODE_BDF(linear_solver = :GMRES)));
        abstol         = 1e-6,
        reltol         = 1e-5,
        save_positions = ((false, false), (false, true)),
        verbose        = false,
        progress,
        plasticity_obs)
    end

    if use_protocol
        @warn (data_protocol[k, :protocol])
    else
        @warn ("Manual Mode Finished")
    end
end

function format_table(dp, selected_idx=1)
        df = copy(dp)
        DataFrames.insertcols!(df, 1, "ID" => 1:DataFrames.nrow(df))

        # Create table header with a sticky position
        headers = [Bonito.DOM.th(string(name), style="position: sticky; top: 0; background-color: #2c3e50; color: white; padding: 8px; border: 1px solid #ddd; z-index: 10;") for name in names(df)]
        thead = Bonito.DOM.thead(Bonito.DOM.tr(headers...))

        # Create table rows with dynamic highlighting for the selected experiment
        rows = map(enumerate(eachrow(df))) do (i, row)
            # Highlight the selected row in light blue, otherwise alternate colors
            bg_color = i == selected_idx ? "#dbeafe" : (iseven(i) ? "#f9f9f9" : "#ffffff")
            font_weight = i == selected_idx ? "bold" : "normal"
            text_color = i == selected_idx ? "#1e40af" : "black"

            cells = [Bonito.DOM.td(string(val), style="padding: 6px 8px; border: 1px solid #ddd; white-space: nowrap;") for val in row]
            Bonito.DOM.tr(cells..., style="background-color: $bg_color; font-weight: $font_weight; color: $text_color; transition: background-color 0.3s;")
        end
        tbody = Bonito.DOM.tbody(rows...)

        # Return the complete HTML table
        return Bonito.DOM.table(thead, tbody, style="width: 100%; border-collapse: collapse; font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; font-size: 12px; text-align: center;")
end

function get_valid_val(tb, min_v, max_v, name)
    val = tryparse(Float64, tb.value[])
    if isnothing(val)
        return nothing, "Impossible: Value of $name is not a valid number."
    elseif val < min_v || val > max_v
        return nothing, "Impossible: Value of $name must be between $min_v and $max_v."
    end
    return val, ""
end

function get_valid_times(tb, name)
    s = strip(tb.value[])
    if isempty(s)
        return Float64[], "" # Autorise une liste vide (ex: aucun Post)
    end
    cleaned = replace(s, "," => " ")
    parts = split(cleaned, keepempty=false)
    times = Float64[]
    for p in parts
        v = tryparse(Float64, p)
        if isnothing(v)
            return Float64[], "Error: Invalid number '$p' in $name."
        end
        push!(times, v)
    end
    return times, ""
end

function build_app(; fontsize = 20)
    app = Bonito.App() do session
        # ── Observables ────────────────────────────────────────────
        p_bar          = Observable{Float64}(0.0)
        plasticity_obs = Observable((
        t      = zeros(1),
        camk2  = zeros(1),
        can    = zeros(1),
        weight = zeros(1)))
        data = Observable(Float64[])

        # ── Widgets Bonito ───────────────
        protocol_options = collect(SynapseElife._protocols)

        # Find the index of TigaretMellor16 in the list
        idx_default = findfirst(==("TigaretMellor16"), protocol_options)
        # Fallback to 1 if not found
        idx_default = isnothing(idx_default) ? 1 : idx_default
        # Force this index on startup
        menu_prot    = Bonito.Dropdown(protocol_options; index=idx_default)

        n_exp_init   = size(dataProtocol(protocol_options[idx_default]), 1)
        menu_exp     = Bonito.Dropdown(collect(1:n_exp_init))
        btn_simulate = Bonito.Button("Run")
        btn_erase    = Bonito.Button("Erase")
        btn_stop        = Bonito.Button("STOP")
        btn_export_csv  = Bonito.Button("Export CSV")
        btn_export_plot = Bonito.Button("Export Plot")

        # ---  MANUAL MODE WIDGETS ---
        mode_manuel = Bonito.Checkbox(false)

        # Text fields (with default values)
        tb_ca    = Bonito.TextField("2.5")
        tb_mg    = Bonito.TextField("1.0")
        tb_temp  = Bonito.TextField("35.0")

        # New exact time textfields (Romain's idea)
        tb_pre_times  = Bonito.TextField("0")
        tb_post_times = Bonito.TextField("10")

        tb_pls   = Bonito.TextField("300")
        tb_freq  = Bonito.TextField("5.0")

        # Dropdowns for Age and Order
        menu_age   = Bonito.Dropdown(["Young", "Adult"]; index=2) # 2 = Adult

        # Error messages
        error_msg = Bonito.Observable("")
        confirm_stop_msg = Bonito.Observable("")
        current_task = Ref{Task}() # Store currently running simulation

        # ── WGLMakie Figure ────────────────────────────────────────
        fig = Figure(;size = (1100, 750), fontsize)

        # Protocol text
        ax12 = Axis(fig[1, 1:3],
        title  = "Experiment information — List of protocols",
        height = 130)
        hidespines!(ax12); hidedecorations!(ax12)
        xlims!(ax12, 0, 15); ylims!(ax12, 0, 2)

        data_protocol_str = Observable(format_table(dataProtocol(protocol_options[idx_default])))

        # Simulation axes
        ax21 = Axis(fig[2, 1],
        xlabel = "time (s)", ylabel = "μM",
        ylabelsize = fontsize / 1.5, xlabelsize = fontsize / 1.5,
        yticklabelsize = fontsize / 2, xticklabelsize = fontsize / 2)
        ylims!(ax21, 0, 34); xlims!(ax21, 0, 200)

        ax22 = Axis(fig[2, 2],
        ylabel = "CaMKII (μM)", xlabel = "CaN (μM)",
        ylabelsize = fontsize / 1.5, xlabelsize = fontsize / 1.5,
        yticklabelsize = fontsize / 2, xticklabelsize = fontsize / 2)

        ax23 = Axis(fig[2, 3],
        ylabel = "weight change (%)", xlabel = "time (s)",
        ylabelsize = fontsize / 1.5, xlabelsize = fontsize / 1.5,
        yticklabelsize = fontsize / 2, xticklabelsize = fontsize / 2)
        xlims!(ax23, 0, 200); ylims!(ax23, -150, 150)

        # Progress bar
        ax_prog = Axis(fig[3, 1:3], title = "Progress", height = 18, yticks = [0])
        pbar_rect = map(p_bar) do x
            Rect(0.0, 0.0, x, 1.0)
        end
        poly!(ax_prog, Rect(0.0, 0.0, 100.0, 1.0),
        color = :white, strokecolor = :red, strokewidth = 1)
        poly!(ax_prog, pbar_rect, color = :red)

        # Live curves
        lines!(ax21,
        lift(x -> x.t ./ 1000, plasticity_obs),
        lift(x -> x.camk2,     plasticity_obs); label = "CaMKII")
        lines!(ax21,
        lift(x -> x.t ./ 1000, plasticity_obs),
        lift(x -> x.can,       plasticity_obs); label = "CaN")
        axislegend(ax21)

        lines!(ax22,
        lift(x -> x.can,   plasticity_obs),
        lift(x -> x.camk2, plasticity_obs))

        lines!(ax23,
        lift(x -> x.t ./ 1000, plasticity_obs),
        lift(x -> x.weight,    plasticity_obs))

        # LTP / LTD regions
        poly!(ax22,
        Point2f[(P[1], P[2]) for P in SynapseParams().LTP_region.vertices];
        color = ("red",  0.1))
        poly!(ax22,
        Point2f[(P[1], P[2]) for P in SynapseParams().LTD_region.vertices];
        color = ("blue", 0.1))

        # ── Callbacks ─────────────────────────────────────────────

        # Update table and reset experiment number when protocol changes
        on(session, menu_prot.value) do s
            dp = dataProtocol(s)
            menu_exp.options[] = collect(1:size(dp, 1))
            menu_exp.value[] = 1
            data_protocol_str[] = format_table(dp, 1)
        end

        # Dynamically highlight the selected experiment row
        on(session, menu_exp.value) do exp_idx
            if !isnothing(exp_idx)
                dp = dataProtocol(menu_prot.value[])
                data_protocol_str[] = format_table(dp, exp_idx)
            end
        end

        on(session, btn_simulate.value) do _
            use_manual = mode_manuel.value[]

            if !use_manual
                error_msg[] = "" # Clear error if in normal mode

                # --- PROTOCOL MODE ---
                k    = parse(Int, string(menu_exp.value[]))
                prot = menu_prot.value[]
                dp   = dataProtocol(prot)
                println("Starting Protocol: ", prot, " (#", k, ")")

                current_task[] = errormonitor(@async begin
                    try
                        run_simulation!(
                        p_bar, data, dp, k;
                        plasticity_obs,
                        use_protocol=true
                        )
                    catch e
                        if e isa InterruptException
                            println("Simulation cancelled by the user.")
                        else
                            rethrow(e)
                        end
                    end
                end)

            else
                # --- MANUAL MODE WITH ALERTS --
                # Check all fields one by one
                v_ca, err = get_valid_val(tb_ca, 2.5, 4.0, "Ca2+")
                if err != ""; error_msg[] = err; return; end
                v_mg, err = get_valid_val(tb_mg, 0.8, 3.0, "Mg2+")
                if err != ""; error_msg[] = err; return; end
                v_temp, err = get_valid_val(tb_temp, 20.0, 40.0, "Temp")
                if err != ""; error_msg[] = err; return; end

                v_pre_times, err = get_valid_times(tb_pre_times, "Pre times")
                if err != ""; error_msg[] = err; return; end
                v_post_times, err = get_valid_times(tb_post_times, "Post times")
                if err != ""; error_msg[] = err; return; end

                v_pls, err = get_valid_val(tb_pls, 1.0, 10000.0, "Pulses")
                if err != ""; error_msg[] = err; return; end
                v_freq, err = get_valid_val(tb_freq, 0.1, 100.0, "Freq")
                if err != ""; error_msg[] = err; return; end

                error_msg[] = ""
                println("▶ Starting Manual Mode...")

                current_task[] = errormonitor(@async begin
                    try
                        run_simulation!(
                            p_bar, data, dataProtocol("TigaretMellor16"), 1;
                            plasticity_obs = plasticity_obs,
                            use_protocol = false,
                            ca_ext = v_ca, mg = v_mg, temp = v_temp,
                            pre_times = v_pre_times, post_times = v_post_times,
                            pulse = Int(v_pls), freq = v_freq, age = menu_age.value[]
                        )
                    catch e
                        if e isa InterruptException
                            println("Simulation cancelled by the user.")
                        else
                            rethrow(e)
                        end
                    end
                end)
            end
        end

        on(session, btn_erase.value) do _
            p_bar[]          = 0.0
            plasticity_obs[] = (t = zeros(1), camk2  = zeros(1),
            can = zeros(1), weight  = zeros(1))
        end

        on(session, btn_export_csv.value) do _
            # Retrieve data
            obs_data = plasticity_obs[]
            # Store in a dataframe
            data_to_export = DataFrame(
            Time_s = obs_data.t ./ 1000, # Convert to seconds
            CaMKII_uM = obs_data.camk2,
            CaN_uM = obs_data.can,
            WeightChange_percent = obs_data.weight
            )

            # Save file
            CSV.write("resultats_synapse.csv", data_to_export)
            println("Data successfully exported to 'resultats_synapse.csv'")
        end

        # We work explicitly with CairoMakie (full namespace)
        # to avoid interfering with active WGLMakie

        on(session, btn_export_plot.value) do _
            # Retrieve current simulation data
            obs = plasticity_obs[]
            t_sec = obs.t ./ 1000

            # 1. Create and save Plot 1: Enzymes Concentration
            fig1 = CairoMakie.Figure(size = (800, 500))
            ax1  = CairoMakie.Axis(fig1[1, 1];
                xlabel = "Time (s)", ylabel = "Concentration (μM)")
            CairoMakie.lines!(ax1, t_sec, obs.camk2; label = "CaMKII", linewidth = 2)
            CairoMakie.lines!(ax1, t_sec, obs.can;   label = "CaN",    linewidth = 2)
            CairoMakie.axislegend(ax1)
            CairoMakie.save("plot_1_enzymes.png", fig1; px_per_unit = 2)

            # 2. Create and save Plot 2: Phase space (CaMKII vs CaN)
            fig2 = CairoMakie.Figure(size = (600, 500))
            ax2  = CairoMakie.Axis(fig2[1, 1];
                xlabel = "CaN (μM)", ylabel = "CaMKII (μM)")
            ltp  = SynapseParams().LTP_region.vertices
            ltd  = SynapseParams().LTD_region.vertices
            CairoMakie.poly!(ax2, CairoMakie.Point2f[(P[1], P[2]) for P in ltp];
               color = (:red,  0.1))
            CairoMakie.poly!(ax2, CairoMakie.Point2f[(P[1], P[2]) for P in ltd];
               color = (:blue, 0.1))
            CairoMakie.lines!(ax2, obs.can, obs.camk2; linewidth = 2)
            CairoMakie.save("plot_2_phase.png", fig2; px_per_unit = 2)

            # 3. Create and save Plot 3: Weight change
            fig3 = CairoMakie.Figure(size = (800, 500))
            ax3  = CairoMakie.Axis(fig3[1, 1];
                xlabel = "Time (s)", ylabel = "Weight change (%)")
            CairoMakie.lines!(ax3, t_sec, obs.weight; linewidth = 2)
            CairoMakie.save("plot_3_weight.png", fig3; px_per_unit = 2)

            println("The 3 plots have been successfully exported as separate files")
        end
        
        on(session, btn_stop.value) do _
            # We stop the simulation if it is still running
            if isassigned(current_task) && !istaskdone(current_task[])
                Base.throwto(current_task[], InterruptException())
                confirm_stop_msg[] = "Simulation stopped successfully."
            else
                confirm_stop_msg[] = "No simulation is currently running."
            end
            
            # The message disappears by itself after 3 seconds
            @async begin
                sleep(3)
                confirm_stop_msg[] = ""
            end
        end
        
        # ── DOM ───────────────────────────────────────────────────

        # Base CSS styles for cards
            card_style = "background: white; padding: 20px; border-radius: 12px; box-shadow: 0 4px 12px rgba(0,0,0,0.08); margin-bottom: 20px;"
            title_style = "margin-top: 0; margin-bottom: 15px; color: #2c3e50; font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; font-size: 18px; border-bottom: 2px solid #ecf0f1; padding-bottom: 8px;"
            
            return Bonito.DOM.div(
                Bonito.DOM.h1("SynapseElife Studio", style="color: #2c3e50; font-family: 'Segoe UI', sans-serif; text-align: center; margin-bottom: 30px; letter-spacing: 1px;"),
                
                # Block 1: Simulation Controls
                Bonito.DOM.div(
                    Bonito.DOM.h3("1. Simulation Controls", style=title_style),
                    Bonito.DOM.div(
                        Bonito.DOM.span("Protocol: ", style="font-weight: bold;"), menu_prot, 
                        Bonito.DOM.span("  Exp #: ", style="font-weight: bold; margin-left: 10px;"), menu_exp, 
                        Bonito.DOM.span(" | ", style="color: #ddd; margin: 0 15px;"),
                        btn_simulate, btn_stop, btn_erase,
                        Bonito.DOM.span(" | ", style="color: #ddd; margin: 0 15px;"), 
                        btn_export_csv, btn_export_plot;
                        style="display: flex; gap: 10px; align-items: center; flex-wrap: wrap;"
                    );
                    style=card_style
                ),

                # Alert Messages
                Bonito.DOM.div(confirm_stop_msg, style="color: #d97706; font-weight: bold; margin-bottom: 5px; text-align: center; font-family: sans-serif;"),
                Bonito.DOM.div(error_msg, style="color: red; font-weight: bold; margin-bottom: 15px; text-align: center; font-family: sans-serif;"),

                # Block 2: Manual Mode
                Bonito.DOM.div(
                    Bonito.DOM.h3("2. Manual Parameters Override", style=title_style),
                    Bonito.DOM.div(
                        Bonito.DOM.span("Enable Manual Mode: ", style="font-weight: bold; color: #c0392b;"), mode_manuel;
                        style="margin-bottom: 20px; background: #fdf2f0; padding: 10px; border-radius: 6px; display: inline-block;"
                    ),
                    Bonito.DOM.div(
                        # Grouping inputs in a CSS Grid for perfect alignment
                        Bonito.DOM.div(Bonito.DOM.span("Ca²⁺ (μM) [2.5-4]:", style="font-size:12px; color:#7f8c8d;"), tb_ca, style="display:flex; flex-direction:column; gap:4px;"),
                        Bonito.DOM.div(Bonito.DOM.span("Mg²⁺ (mM) [0.8-3]:", style="font-size:12px; color:#7f8c8d;"), tb_mg, style="display:flex; flex-direction:column; gap:4px;"),
                        Bonito.DOM.div(Bonito.DOM.span("Temp (°C):", style="font-size:12px; color:#7f8c8d;"), tb_temp, style="display:flex; flex-direction:column; gap:4px;"),
                        Bonito.DOM.div(Bonito.DOM.span("Pre times (ms):", style="font-size:12px; color:#7f8c8d; font-weight:bold;"), tb_pre_times, style="display:flex; flex-direction:column; gap:4px;"),
                        Bonito.DOM.div(Bonito.DOM.span("Post times (ms):", style="font-size:12px; color:#7f8c8d; font-weight:bold;"), tb_post_times, style="display:flex; flex-direction:column; gap:4px;"),
                        Bonito.DOM.div(Bonito.DOM.span("Repetitions (Pulses):", style="font-size:12px; color:#7f8c8d;"), tb_pls, style="display:flex; flex-direction:column; gap:4px;"),
                        Bonito.DOM.div(Bonito.DOM.span("Frequency (Hz):", style="font-size:12px; color:#7f8c8d;"), tb_freq, style="display:flex; flex-direction:column; gap:4px;"),
                        Bonito.DOM.div(Bonito.DOM.span("Subject Age:", style="font-size:12px; color:#7f8c8d;"), menu_age, style="display:flex; flex-direction:column; gap:4px;");
                        style="display: grid; grid-template-columns: repeat(auto-fit, minmax(130px, 1fr)); gap: 20px;"
                    );
                    style=card_style
                ),

                # Block 3: Protocol Table
                Bonito.DOM.div(
                    Bonito.DOM.h3("Protocol Details (Selected row is highlighted)", style=title_style),
                    Bonito.DOM.div(data_protocol_str, style="overflow-y: auto; overflow-x: auto; max-height: 250px; border: 1px solid #ecf0f1; border-radius: 8px;"),
                    style=card_style
                ),

                # Block 4: Simulation Plots
                Bonito.DOM.div(
                    fig,
                    style="background: white; padding: 20px; border-radius: 12px; box-shadow: 0 4px 12px rgba(0,0,0,0.08); display: flex; justify-content: center;"
                );
                
                # Global page style
                style="background-color: #f4f7f6; padding: 40px; min-height: 100vh; box-sizing: border-box;"
            )

            
            # The Figure
            fig
    end     
    return app
end


# ── Launch ─────────────────────────────────────────────────────
app    = build_app()
server = Bonito.Server(app, "127.0.0.1", 8081)

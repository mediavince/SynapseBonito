# Synapse E-Life Studio

# SynapseBonito

SynapseElife Studio — Interactive GUI  (v6 — tabbed layout)

Built with Bonito (web UI) + WGLMakie (live plots) + CairoMakie (export)

##  ARCHITECTURE OVERVIEW

    1. Helper functions   — input parsing, event construction
    2. Simulation engine  — run_one_simulation! and run_simulation!
    3. GUI helpers        — HTML table formatter
    4. build_app()        — the full Bonito application (widgets, plots, callbacks)

Bonito registers all widgets by their DOM id at page load time.

Hiding widgets with display:none prevents Bonito from finding them → crashes.

Solution: ALL widgets are rendered into the DOM once at startup (inside invisible wrapper divs).  Tab switching is handled by a Julia-side Observable{String} called `active_tab`.  The visible tab content is produced  by lift(active_tab) which returns a DOM node containing only the active panel's content.  Inactive panels are rendered as empty divs.

##  SIMULATION MODES

    • Protocol : reads parameters from a pre-defined experiment table
    • STDP     : loops over Δt values, builds the full STDP curve
    • HFS/LFS  : pre-synaptic train at a user-defined frequency
    • TBS      : theta-burst stimulation (N bursts × P pulses @ intra Hz)
    • Custom   : fully custom pre/post spike timestamps

##  EXPORT

WGLMakie renders in the browser via WebGL — it cannot write PNG files.

CairoMakie (server-side) is used exclusively for file export.

Tabular data and simulation results are exported as timestamped .xlsx files.

## DEV

# Pre-requisite

```shell
winget install --name Julia --id 9NJNWW8PVKMN -e -s msstore
```

# Launch Program locally

```shell
julia --project=@. app_bonito_vf.jl
```
 or
```shell
julia --project=@synapsebonito app_bonito_vf.jl
```

## NB

juliaup is the official recommended method — it's the only approach the Julia docs endorse, and it avoids all the Chocolatey problems (PATH conflicts, stale versions, no easy switching). 

Install (one-liner):

winget install --name Julia --id 9NJNWW8PVKMN -e -s msstore

This installs juliaup (the version manager) from the Microsoft Store, which automatically pulls the latest stable Julia.  No PATH edits, no Chocolatey, no conflicts.

Why it's better than the alternatives:

Method	Problem
Chocolatey (choco install julia)	Puts julia.exe in C:\ProgramData\chocolatey\bin\ — conflicts with any other Julia on PATH, pins to whatever version choco last updated, no easy version switching
Manual installer (.exe from julialang.org)	Works, but you manually re-download + re-install for every upgrade, and the PATH entry is tied to a specific version folder
juliaup ✅	Single shim on PATH (WindowsApps\julia.exe), manages multiple versions, juliaup update upgrades in-place, juliaup add 1.11 / juliaup default 1.11 for pinning

After installing, manage versions with:

juliaup list              # show installed versions
juliaup add 1.11          # install a specific version
juliaup default 1.11      # set as default
juliaup update            # update to latest patch

Key rule: never have more than one Julia source on your PATH. If you had Chocolatey's Julia before, choco uninstall julia -y first (as you already did), then install juliaup. The juliaup shim lives in WindowsApps and won't conflict with anything else.

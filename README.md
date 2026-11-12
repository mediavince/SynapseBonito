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

## CI/CD

use the formatter before pushing

For formatting, don't do it by hand — use JuliaFormatter:

julia -e 'using Pkg; Pkg.add("JuliaFormatter")'
julia -e 'using JuliaFormatter; format_file("app_bonito_vf.jl")'

Or for the whole project:

julia -e 'using JuliaFormatter; format_dir(".")'

This will fix all the inconsistent indentation (your script mixes 4-space, 8-space, and tab indentation in places). Add a .JuliaFormatter.toml to your project root to lock the style:

indent = 4
margin = 120

Then format_dir("src") will produce consistent output.


## Run with Docker Compose

Requirements: Docker Engine with Docker Compose v2, or Docker Desktop.

Edit .env-dist with your preferred variables and save as .env

Run these commands from the repository root:

```shell
docker compose up -d --build
docker compose logs -f --tail 50 synapse_elife_studio
```

Wait for `Server started`, then open http://localhost:8081.

The image uses Julia 1.11.9. Dependencies are installed from
`src/Project.toml` and `src/Manifest.toml` and precompiled during the build.
The first build may take several minutes.

Press Ctrl+C to stop following logs; the application keeps running.

### Stop and restart

```shell
docker compose stop
docker compose start
```

To rebuild after changing application files:

```shell
docker compose up -d --build
```

Application code is included in the image, so source changes require a rebuild.

### Exported files

The export buttons write Excel and PNG files to `/app/exports`.
This directory uses the persistent Docker volume `synapse_exports`.

Files are not automatically downloaded by the browser.

To copy exports to the host:

```shell
mkdir -p exports-recuperes
docker compose cp synapse_elife_studio:/app/exports/. ./src/exports/
```

The volume survives container replacement and `docker compose down`.
Do not use `docker compose down -v` unless you intend to delete the exports.

## Remote deployment

The default Compose configuration exposes the application only on the
host's loopback address, `127.0.0.1:8081`.

For remote access, configure an HTTPS reverse proxy and set `PUBLIC_URL`
to the address used by visitors. For example, create a `.env` file in the
repository root:

```dotenv
# Example environment variables
JULIA_ENV=production

# Domain like domain.lan or localhost
DOMAIN=synapse.example.org

# NGINX Listen Port (exposed to host)
NGINX_PORT=8081

# Internal port for Julia app
JULIA_HOST=0.0.0.0
JULIA_PORT=8000

# SSL Certificates (paths inside the container) uncomment if needed
#SSL_CERT=/etc/nginx/ssl/cert.pem
#SSL_KEY=/etc/nginx/ssl/key.pem
```

Then apply the configuration:

The dockerfile already prepares the required setup with the nginx.conf template that has placeholders to input .env related vars.

```shell
docker compose up -d --build
```

The reverse proxy must forward HTTP requests, static assets and WebSocket
connections to the application. A proxy running directly on the same host
can forward to `http://127.0.0.1:8081`.

`JULIA_HOST` and `JULIA_PORT` control the listener inside the container.
`DOMAIN` controls the external address advertised by Bonito.

Using `localhost` as the public URL will not work for remote visitors.

Build and test the image for the target server's CPU architecture.
Successful testing on a local machine does not validate a remote deployment.

## Run locally without Docker

With Juliaup installed, install the matching Julia version:

```shell
juliaup add 1.11.9
```

From the repository root, install dependencies once and launch:

```shell
cd src
julia +1.11.9 --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'
julia +1.11.9 --project=. app_bonito_vf.jl
```

Open http://localhost:8081. Stop any Docker instance using this port first.

Repeat dependency installation after changes to Project.toml or Manifest.toml.

## Validation and CI

The GitHub Actions workflow checks formatting, builds the Docker image
on Linux and checks that the application returns a non-empty HTTP page.
Application logs are collected, including when the startup check fails.

This HTTP check does not validate browser interactions or simulations.

Before deploying, manually check:

- Graphs appear and controls respond.
- A simulation completes and displays its results.
- Excel and PNG exports are generated and can be retrieved.
- The application works after a restart.
- Exports survive container replacement.

Multi-user capacity, resource limits, access control and the target
server's HTTPS/WebSocket configuration require separate validation.

## Nginx build with ssl certs

Ensure the .env (edited from .env-dist) contains ssl vars for https domain; uncomment ssl related lines in docker-compose.yml

FROM julia:1.11.9

WORKDIR /app

RUN apt-get update && apt-get install -y --no-install-recommends \
    libgl1 \
    libglu1-mesa \
    libx11-6 \
    libxt6 \
    libxrender1 \
    libxext6 \
    xvfb \
    xauth \
    libcairo2 \
    libpangocairo-1.0-0 \
    fontconfig \
    fonts-dejavu-core \
    && rm -rf /var/lib/apt/lists/*

ENV GKS_WSTYPE=png \
    JULIA_PKG_PRECOMPILE_AUTO=0 \
    JULIA_NUM_PRECOMPILE_TASKS=2

# Installer les versions enregistrees dans le Manifest.
COPY src/Project.toml src/Manifest.toml ./

RUN julia --project=/app -e 'using Pkg; Pkg.instantiate()'

# Precompiler pendant la construction, avec un ecran virtuel disponible.
RUN xvfb-run -a -e /dev/stderr \
    --server-args="-screen 0 1920x1080x24" \
    julia --project=/app -e 'using Pkg; Pkg.precompile()'

COPY src/ ./

EXPOSE 8081

CMD ["xvfb-run", "-a", "-e", "/dev/stderr", "--server-args=-screen 0 1920x1080x24", "julia", "--project=/app", "app_bonito_vf.jl"]
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
    nginx \
    && rm -rf /var/lib/apt/lists/*

ENV GKS_WSTYPE=png \
    JULIA_PKG_PRECOMPILE_AUTO=0 \
    JULIA_NUM_PRECOMPILE_TASKS=2

COPY src/Project.toml ./

# RUN julia -e 'using Pkg; Pkg.add(["Bonito","CSV","CairoMakie","DataFrames","Dates","Hyperscript","JuliaFormatter","PiecewiseDeterministicMarkovProcesses","PrettyTables","ProgressMeter","Random","Revise","Sundials","Symbolics","SynapseElife","WGLMakie","XLSX"]);'

# Remove default NGINX config to avoid conflicts
RUN rm -f /etc/nginx/conf.d/default.conf

# Copy nginx.conf to the container
COPY nginx.conf /etc/nginx/conf.d/default.conf

ARG DOMAIN=localhost
ARG NGINX_PORT=8081
ARG JULIA_PORT=8000

RUN sed -i \
    -e "s/__DOMAIN__/${DOMAIN}/g" \
    -e "s/__NGINX_PORT__/${NGINX_PORT}/g" \
    -e "s/__JULIA_PORT__/${JULIA_PORT}/g" \
    /etc/nginx/conf.d/default.conf && \
    nginx -t

RUN julia --project=/app -e 'using Pkg; Pkg.instantiate()'

RUN xvfb-run -a -e /dev/stderr \
    --server-args="-screen 0 1920x1080x24" \
    julia --project=/app -e 'using Pkg; Pkg.precompile()'

COPY src/ ./

EXPOSE 8081

CMD ["sh", "-c", "nginx -g 'daemon off;' & exec xvfb-run -a -e /dev/stderr --server-args='-screen 0 1920x1080x24' julia --project=/app app_bonito_vf.jl"]

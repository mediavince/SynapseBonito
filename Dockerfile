FROM julia:1.11

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

ENV GKS_WSTYPE=png

RUN julia -e 'using Pkg; Pkg.add(["Bonito","CSV","CairoMakie","DataFrames","DotEnv","GLMakie","Hyperscript","PiecewiseDeterministicMarkovProcesses","PrettyTables","ProgressMeter","Revise","Sundials","Symbolics","SynapseElife","WGLMakie","XLSX"]);'

# Instantiate, will instantiate properly upon first call, for now as build fails
#RUN xvfb-run --server-args="-screen 0 1920x1080x24" julia -e 'using Pkg; Pkg.precompile(); Pkg.gc()'

COPY ./src /app/

CMD ["xvfb-run", "--server-args=-screen 0 1920x1080x24", "julia", "--project=@.", "app_bonito_vf.jl"]

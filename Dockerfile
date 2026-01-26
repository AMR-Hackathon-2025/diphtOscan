FROM condaforge/mambaforge:latest

LABEL maintainer="Martin RETHORET-PASTY <martin.rethoret-pasty@pasteur.fr>"
LABEL description="diphtOscan: A tool for characterising virulence and resistance in Corynebacterium diphtheriae"
LABEL org.opencontainers.image.source="https://github.com/AMR-Hackathon-2025/diphtOscan"
LABEL org.opencontainers.image.licenses="GPL-3.0-or-later"

SHELL ["/bin/bash", "-o", "pipefail", "-c"]

WORKDIR /opt/diphtoscan

# Copy environment file first for layer caching
COPY environment.yml .

# Create conda environment
RUN mamba env create -f environment.yml && \
    mamba clean -afy && \
    rm -rf /opt/conda/pkgs/*

# Copy source code
COPY pyproject.toml README.md LICENSE ./
COPY diphtoscan/ ./diphtoscan/

# Install package (set version for setuptools-scm since .git is not copied)
ARG VERSION=1.7.0
RUN SETUPTOOLS_SCM_PRETEND_VERSION=${VERSION} conda run -n diphtoscan pip install . --no-deps --no-cache-dir

# Activate environment
ENV PATH="/opt/conda/envs/diphtoscan/bin:${PATH}"
ENV CONDA_DEFAULT_ENV="diphtoscan"

# Optional: Update database at build time (uncomment to include)
# RUN conda run -n diphtoscan diphtoscan -u

# Create non-root user
RUN useradd --create-home --shell /bin/bash --uid 1000 diphtoscan && \
    mkdir -p /data /output && \
    chown -R diphtoscan:diphtoscan /data /output /opt/diphtoscan

WORKDIR /data
USER diphtoscan

ENTRYPOINT ["diphtoscan"]
CMD ["--help"]

HEALTHCHECK --interval=30s --timeout=10s --start-period=5s --retries=3 \
    CMD diphtoscan --version || exit 1

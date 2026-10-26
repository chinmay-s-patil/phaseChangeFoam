# syntax=docker/dockerfile:1
ARG OPENFOAM_VERSION=2412
FROM opencfd/openfoam-default:${OPENFOAM_VERSION}

USER root

# Install Python 3 and verification plotting dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    python3 \
    python3-numpy \
    python3-matplotlib \
    python3-pil \
    && rm -rf /var/lib/apt/lists/*

# Copy phaseChangeFoam source into container
ENV FOAM_PROJECT_DIR=/usr/local/src/phaseChangeFoam
WORKDIR ${FOAM_PROJECT_DIR}
COPY . ${FOAM_PROJECT_DIR}

# Compile solvers and utilities into OpenFOAM environment and system path
RUN bash -c '\
    FOAM_BASHRC="/usr/lib/openfoam/openfoam'${OPENFOAM_VERSION}'/etc/bashrc"; \
    if [ ! -f "$FOAM_BASHRC" ]; then \
        FOAM_BASHRC=$(find /usr/lib/openfoam /opt -name bashrc 2>/dev/null | grep etc/bashrc | head -n 1); \
    fi; \
    echo "Using OpenFOAM bashrc: $FOAM_BASHRC"; \
    source "$FOAM_BASHRC" && \
    ./Allwmake && \
    cp -f ${FOAM_USER_APPBIN}/* ${FOAM_APPBIN}/ 2>/dev/null || true && \
    cp -f ${FOAM_USER_APPBIN}/* /usr/local/bin/ 2>/dev/null || true && \
    chmod 755 /usr/local/bin/phaseChangeMultiRegion* /usr/local/bin/writePCDict 2>/dev/null || true'

# Smoke test binaries to guarantee successful compilation and linking
RUN bash -c '\
    FOAM_BASHRC="/usr/lib/openfoam/openfoam'${OPENFOAM_VERSION}'/etc/bashrc"; \
    if [ ! -f "$FOAM_BASHRC" ]; then \
        FOAM_BASHRC=$(find /usr/lib/openfoam /opt -name bashrc 2>/dev/null | grep etc/bashrc | head -n 1); \
    fi; \
    source "$FOAM_BASHRC" && \
    phaseChangeMultiRegionFoam -help && \
    phaseChangeMultiRegionSimpleFoam -help && \
    writePCDict -help'

# Switch to default non-root openfoam user
USER openfoam
WORKDIR /home/openfoam

CMD ["/bin/bash"]

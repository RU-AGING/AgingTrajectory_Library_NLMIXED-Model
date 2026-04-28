# syntax=docker/dockerfile:1.6
#
# Traj2 — GBTM Macro Library
# https://github.com/RU-AGING/AgingTrajectory_Library_NLMIXED-Model
#
# This image contains the Traj2 SAS macros, simulators, plotting utilities,
# and entrypoint scripts. It does NOT contain SAS itself : SAS is proprietary
# and cannot be redistributed. Users must mount their licensed SAS 9.4+
# installation read-only at /opt/sas at runtime.

FROM ubuntu:22.04

LABEL org.opencontainers.image.title="Traj2"
LABEL org.opencontainers.image.description="Group-Based Trajectory Modeling SAS macro library (Ordinal-Probit, truncated-normal continuous, ZIP, ZINB)"
LABEL org.opencontainers.image.source="https://github.com/RU-AGING/AgingTrajectory_Library_NLMIXED-Model"
LABEL org.opencontainers.image.licenses="CC-BY-4.0"
LABEL org.opencontainers.image.vendor="CHAO Lab, Rutgers University"

ENV DEBIAN_FRONTEND=noninteractive \
    SAS_HOME=/opt/sas \
    GBTM_HOME=/opt/gbtm \
    GBTM_MACROS=/opt/gbtm/macros \
    GBTM_DATA=/opt/gbtm/data \
    GBTM_OUTPUT=/opt/gbtm/output

# Minimal OS deps. SAS itself provides its own runtime; we only need
# bash, file utilities, and a small text-processing toolchain for the
# entrypoint script.
RUN apt-get update && apt-get install -y --no-install-recommends \
        bash \
        ca-certificates \
        coreutils \
        findutils \
        gawk \
        grep \
        sed \
        tzdata \
    && rm -rf /var/lib/apt/lists/*

# Layout:
#   /opt/gbtm/macros   — the SAS source files copied from the repo
#   /opt/gbtm/data     — user mounts their input data here
#   /opt/gbtm/output   — user mounts an output directory here
#   /opt/sas           — user mounts their licensed SAS install here (RO)
RUN mkdir -p ${GBTM_MACROS} ${GBTM_DATA} ${GBTM_OUTPUT}

# Copy the SAS source files. We copy the whole repo's docx/rtf/etc. so the
# image stays self-contained even if files are renamed; the entrypoint
# resolves them by name at run time.
COPY discrete_ordinal_final.docx       ${GBTM_MACROS}/
COPY Continuous_2_outcomes.docx        ${GBTM_MACROS}/
COPY poisson_model_complete_35_.rtf    ${GBTM_MACROS}/
COPY ZINB_LatentClass_Trajectories.docx ${GBTM_MACROS}/

# Convenience: extract plain-text .sas wrappers from the .docx/.rtf at build
# time would require pandoc/unoconv. We instead ship the source files as-is
# and document that users include them via SAS file references that handle
# RTF/DOCX. If you have plain .sas exports, drop them into ./macros_sas/
# in the repo and they will be copied below.
COPY macros_sas/ ${GBTM_MACROS}/sas/

# Entrypoint
COPY docker/entrypoint.sh /usr/local/bin/entrypoint.sh
RUN chmod +x /usr/local/bin/entrypoint.sh

WORKDIR ${GBTM_HOME}
ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]
CMD ["help"]

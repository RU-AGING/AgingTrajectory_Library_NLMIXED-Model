# syntax=docker/dockerfile:1.6
#
# Traj2 — GBTM Macro Library
# https://github.com/RU-AGING/AgingTrajectory_Library_NLMIXED-Model
#
# This image contains the Traj2 SAS macros, simulators, plotting utilities,
# and entrypoint scripts. It does NOT contain SAS itself — SAS is proprietary
# and cannot be redistributed. Users must mount their licensed SAS 9.4+
# installation read-only at /opt/sas at runtime.
FROM ubuntu:22.04
LABEL org.opencontainers.image.title="Traj2"
LABEL org.opencontainers.image.description="Group-Based Trajectory Modeling SAS macro library — ordinal-probit and censored-normal continuous (released); ZIP and ZINB (prototypes, under QA)"
LABEL org.opencontainers.image.source="https://github.com/RU-AGING/AgingTrajectory_Library_NLMIXED-Model"
LABEL org.opencontainers.image.licenses="MIT"
LABEL org.opencontainers.image.vendor="CHAO Lab, Rutgers University"
ENV DEBIAN_FRONTEND=noninteractive \
    SAS_HOME=/opt/sas \
    GBTM_HOME=/opt/gbtm \
    GBTM_MACROS=/opt/gbtm/macros \
    GBTM_DATA=/opt/gbtm/data \
    GBTM_OUTPUT=/opt/gbtm/output
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
RUN mkdir -p ${GBTM_MACROS}/docs ${GBTM_DATA} ${GBTM_OUTPUT}
COPY traj2_macro.sas             ${GBTM_MACROS}/
COPY traj2_main_macro.sas        ${GBTM_MACROS}/
COPY discrete_ordinal_final.sas  ${GBTM_MACROS}/
COPY continuous_2_outcomes.sas   ${GBTM_MACROS}/
COPY zip_poisson.sas             ${GBTM_MACROS}/
COPY zinb.sas                    ${GBTM_MACROS}/
COPY discrete_ordinal_final.docx        ${GBTM_MACROS}/docs/
COPY Continuous_2_outcomes.docx         ${GBTM_MACROS}/docs/
COPY poisson_model_complete_35_.rtf     ${GBTM_MACROS}/docs/
COPY ZINB_LatentClass_Trajectories.docx ${GBTM_MACROS}/docs/
COPY docker/entrypoint.sh /usr/local/bin/entrypoint.sh
RUN chmod +x /usr/local/bin/entrypoint.sh
WORKDIR ${GBTM_HOME}
ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]
CMD ["help"]

# Traj2: Group-Based Trajectory Modeling SAS Macro Library

[![Docker Image](https://img.shields.io/badge/docker-ghcr.io%2Fru--aging%2Ftraj2-blue?logo=docker)](https://github.com/RU-AGING/AgingTrajectory_Library_NLMIXED-Model/pkgs/container/traj2)
[![Documentation](https://img.shields.io/badge/docs-portal-teal)](https://ru-aging.github.io/AgingTrajectory_Library_NLMIXED-Model/)
[![License: CC BY 4.0](https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)

A native SAS macro library for **Group-Based Trajectory Modeling (GBTM)** built on `PROC NLMIXED`. Developed at the **Community Health and Aging Outcome (CHAO) Lab**, Rutgers University.

Traj2 fills a gap left by the legacy `PROC TRAJ` for environments such as the CMS Virtual Research Data Center (VRDC) where third-party procedures cannot be installed, and provides outcome distributions and joint-trajectory functionality not available in `PROC TRAJ`.

## Models

| Model | Outcome type | Use cases |
|---|---|---|
| **Ordinal-Probit** | Bounded ordered categorical | ADL/IADL, frailty grades, cognitive screens, pain scales |
| **Truncated-Normal Continuous** (single + joint two-outcome) | Continuous on bounded support | Utilization with per-period caps, biomarkers with floor/ceiling effects |
| **ZIP** *(extension)* | Integer counts with excess zeros | ED visits, hospitalizations, prescription fills |
| **ZINB** *(extension)* | Integer counts with overdispersion | Inpatient admissions, post-acute care days |

## Documentation

Full documentation, model equations, quick-start examples, and download links are in the documentation portal:

**https://ru-aging.github.io/AgingTrajectory_Library_NLMIXED-Model/**

## Quick start (Docker)

The easiest way to use Traj2 is via the prebuilt Docker image. The image bundles all macros, simulators, and entrypoint scripts; it does **not** bundle SAS, which is proprietary — mount your licensed SAS installation at `/opt/sas`:

```bash

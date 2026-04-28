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

The easiest way to use Traj2 is via the prebuilt Docker image. The image bundles all macros, simulators, and entrypoint scripts; it does **not** bundle SAS, which is proprietary, mount your licensed SAS installation at `/opt/sas`:

```bash
docker pull ghcr.io/ru-aging/traj2:1.0.0

docker run --rm \
  -v /path/to/your/sas:/opt/sas:ro \
  -v "$PWD/data":/opt/gbtm/data \
  -v "$PWD/output":/opt/gbtm/output \
  ghcr.io/ru-aging/traj2:1.0.0 \
  run --model=ordinal --nclass=4 --order=2
```

Available models: `ordinal`, `continuous`, `zip`, `zinb`.

## Quick start (SAS, without Docker)

Open one of the `.sas` files in your SAS session and follow the example block at the bottom of each:

- `ordinal_probit.sas` — Ordinal-Probit GBTM
- `continuous_2outcomes.sas` — Joint truncated-normal GBTM for two continuous outcomes
- `zip_poisson.sas` — Zero-Inflated Poisson trajectories
- `zinb.sas` — Zero-Inflated Negative Binomial trajectories

Each file is self-contained, includes an optional simulator, and follows a `set globals → simulate or load data → fit → plot` workflow.

## Repository contents

| File | Purpose |
|---|---|
| `ordinal_probit.sas`, `continuous_2outcomes.sas`, `zip_poisson.sas`, `zinb.sas` | Plain-text SAS macro files (one per model family) |
| `traj2_macro.sas`, `traj2_main_macro.sas` | Earlier macro library files |
| `discrete_ordinal_final.docx`, `Continuous_2_outcomes.docx`, `poisson_model_complete_35_.rtf`, `ZINB_LatentClass_Trajectories.docx` | Documented model code (human-readable reference) |
| `discrete_ordinal_user_guide.docx`, `POission_Final_User_Guide_79__70_.docx` | User guides |
| `TRAJ2_QA_47_.xlsx` | 47 QA test cases |
| `GBTM_Macros_Testing_Brief.pptx` | Analyst testing brief |
| `Dockerfile`, `docker-compose.yml`, `docker/entrypoint.sh` | Docker image build files |
| `BUILD_AND_PUSH.md` | Instructions for rebuilding and publishing the Docker image |
| `index.html` | GitHub Pages documentation portal source |

## Citation

If you use Traj2 in your research, please cite:

> Lin, H., Zafar, A., Xia, W., Jones, B., & Jarrín, O. F. (2026). *Traj2: A Native Macro for Single and Multi-Outcome Group-Based Trajectory Modeling in SAS.* Journal of Statistical Software (in preparation).

## Authors

- **Haiqun Lin** — Methodology, ordinal-probit and truncated-normal models
- **Anum Zafar** — ZIP, ZINB, ordinal-probit; macro library architecture; documentation
- **Weiyi Xia** — Two-outcome continuous model
- **Bobby Jones** — Carnegie Mellon University; methodology consultation
- **Olga F. Jarrín** — Principal Investigator, CHAO Lab

## License

The macro library is released under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). You are free to share and adapt the work with appropriate attribution.

SAS itself is proprietary and is **not** included in this repository or in the Docker image. Users must hold a valid SAS 9.4+ license to run the macros.

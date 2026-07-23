# Traj2: Group-Based Trajectory Modeling SAS Macro Library

[![Docker Image](https://img.shields.io/badge/docker-chao--lab%2Fgbtm--macros-blue?logo=docker)](https://hub.docker.com/r/chao-lab/gbtm-macros)
[![Documentation](https://img.shields.io/badge/docs-portal-teal)](https://ru-aging.github.io/AgingTrajectory_Library_NLMIXED-Model/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A native SAS macro library for **Group-Based Trajectory Modeling (GBTM)** built on `PROC NLMIXED`. Developed at the **Community Health and Aging Outcomes (CHAO) Lab**, Rutgers University.

Traj2 fills a gap left by the legacy `PROC TRAJ` for environments such as the CMS Virtual Research Data Center (VRDC) where externally compiled procedures cannot be installed, and provides outcome distributions and joint-trajectory functionality within an editable, pure-SAS implementation.

## Models

| Model | Outcome type | Use cases |
|---|---|---|
| **Ordinal-Probit** *(released)* | Bounded ordered categorical | ADL/IADL, frailty grades, cognitive screens, pain scales |
| **Censored-Normal Continuous** *(released; single + joint two-outcome)* | Continuous on bounded support | Utilization with per-period caps, biomarkers with floor/ceiling effects |
| **ZIP** *(prototype · under QA · not in the current release)* | Integer counts with excess zeros | ED visits, hospitalizations, prescription fills |
| **ZINB** *(prototype · under QA · not in the current release)* | Integer counts with overdispersion | Inpatient admissions, post-acute care days |

## Documentation

Full documentation, model equations, quick-start examples, and download links are in the documentation portal:

**https://ru-aging.github.io/AgingTrajectory_Library_NLMIXED-Model/**

## Quick start (Docker)

The easiest way to use Traj2 is via the prebuilt Docker image. The image bundles all macros, simulators, and entrypoint scripts; it does **not** bundle SAS, which is proprietary : mount your licensed SAS installation at `/opt/sas`:

```bash
docker pull chao-lab/gbtm-macros:1.0.0

docker run --rm \
  -v /path/to/your/sas:/opt/sas:ro \
  -v "$PWD/data":/opt/gbtm/data \
  -v "$PWD/output":/opt/gbtm/output \
  chao-lab/gbtm-macros:1.0.0 \
  run --model=ordinal --nclass=4 --order=2
```

Released models: `ordinal`, `continuous`. Prototypes (QA testing only): `zip`, `zinb`.

## Quick start (SAS, without Docker)

Open one of the `.sas` files in your SAS session and follow the example block at the bottom of each:

- `ordinal_probit.sas` : Ordinal-Probit GBTM
- `continuous_2outcomes.sas` : Censored-normal GBTM (single + joint two-outcome)
- `zip_poisson.sas` : Zero-Inflated Poisson trajectories (prototype, under QA)
- `zinb.sas` : Zero-Inflated Negative Binomial trajectories (prototype, under QA)

Each file is self-contained, includes an optional simulator, and follows a `set globals → simulate or load data → fit → plot` workflow.

## Repository contents

| File | Purpose |
|---|---|
| `ordinal_probit.sas`, `continuous_2outcomes.sas`, `zip_poisson.sas`, `zinb.sas` | Plain-text SAS macro files (one per model family) |
| `Data_Dictionary_continuous_2outcomes.pdf`, `Data_Dictionary_ordinal_probit.pdf`, `Data_Dictionary_zinb.pdf`, `Data_Dictionary_zip_poisson.pdf` | Data Dictionary |
| `continuous_2outcomes_user_guide.pdf`,`ordinal_probit_user_guide.pdf`,`zinb_user_guide.pdf`, `zinb_user_guide.pdf` | User guides |
| `TRAJ2_QA_test_cases.xlsx` | QA test cases |
| `Dockerfile`, `docker-compose.yml`, `docker/entrypoint.sh` | Docker image build files |
| `BUILD_AND_PUSH.md` | Instructions for rebuilding and publishing the Docker image |
| `index.html` | GitHub Pages documentation portal source |

## Citation

If you use Traj2 in your research, please cite:

> Zafar, A., Xia, W., Lin, H., & Jarrín, O. F. (2026). *Traj2: A Native Macro Library for Single and Multi-Outcome Group-Based Trajectory Modeling in SAS.* Journal of Statistical Software (in preparation).

## Authors

- **Anum Zafar** : Macro library architecture; ordinal-probit; ZIP/ZINB prototypes; documentation
- **Weiyi Xia** : Two-outcome continuous model
- **Haiqun Lin** : Methodology; ordinal-probit and censored-normal models
- **Olga F. Jarrín** : Principal Investigator, CHAO Lab

## License

The macro library **source code** is released under the [MIT License](https://opensource.org/licenses/MIT). The **documentation** (documentation portal and user guides) is released under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).

SAS itself is proprietary and is **not** included in this repository or in the Docker image. Users must hold a valid SAS 9.4+ license to run the macros.

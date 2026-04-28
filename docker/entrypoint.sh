#!/usr/bin/env bash
#
# Traj2 entrypoint
# Usage:
#   docker run … ghcr.io/ru-aging/traj2:VERSION <command> [args]
#
# Commands:
#   help        Show this help (default)
#   validate    Confirm SAS is mounted and runnable
#   list        List available macros
#   run         Run a model. Requires --model and --nclass.
#   shell       Drop into a bash shell inside the container

set -euo pipefail

SAS_BIN=""
SAS_HOME="${SAS_HOME:-/opt/sas}"
GBTM_HOME="${GBTM_HOME:-/opt/gbtm}"
GBTM_MACROS="${GBTM_MACROS:-/opt/gbtm/macros}"
GBTM_DATA="${GBTM_DATA:-/opt/gbtm/data}"
GBTM_OUTPUT="${GBTM_OUTPUT:-/opt/gbtm/output}"

resolve_sas() {
  # Look for the SAS launcher. Common locations:
  #   /opt/sas/SASFoundation/9.4/sas
  #   /opt/sas/sas
  #   $SAS_HOME/sas
  for candidate in \
      "${SAS_HOME}/SASFoundation/9.4/sas" \
      "${SAS_HOME}/sas" \
      "${SAS_HOME}/SASHome/SASFoundation/9.4/sas" \
      "$(command -v sas 2>/dev/null || true)"; do
    if [ -n "${candidate}" ] && [ -x "${candidate}" ]; then
      SAS_BIN="${candidate}"
      return 0
    fi
  done
  return 1
}

cmd_help() {
  cat <<'EOF'
Traj2 — GBTM Macro Library Docker image
========================================

This image bundles the Traj2 SAS macros. SAS itself is NOT included.
You must mount your licensed SAS 9.4+ installation read-only at /opt/sas.

Quick start:

  docker run --rm \
    -v /path/to/your/sas:/opt/sas:ro \
    -v "$PWD/data":/opt/gbtm/data \
    -v "$PWD/output":/opt/gbtm/output \
    ghcr.io/ru-aging/traj2:1.0.0 validate

Commands:
  help                    Show this help
  validate                Confirm SAS is found and runnable
  list                    List available macros
  run --model MODEL       Run a model. MODEL is one of:
                            ordinal | continuous | zip | zinb
      --nclass N          Number of latent classes (required for run)
      --order N           Polynomial order (default 2)
      --T N               Number of time points (default 12)
  shell                   Open a bash shell inside the container

Mounts:
  /opt/sas         (RO)   Your licensed SAS installation
  /opt/gbtm/data          Your input data (BASE_FILE_SRS, etc.)
  /opt/gbtm/output        Where outputs are written

Example — fit a 4-class ordinal model:

  docker run --rm \
    -v /opt/sas:/opt/sas:ro \
    -v "$PWD/data":/opt/gbtm/data \
    -v "$PWD/output":/opt/gbtm/output \
    ghcr.io/ru-aging/traj2:1.0.0 \
    run --model=ordinal --nclass=4 --order=2

Citation:
  Lin, H., Zafar, A., Xia, W., Jones, B., & Jarrín, O. F. (2026).
  Traj2: A Native Macro for Single- and Multi-Outcome Group-Based
  Trajectory Modeling in SAS. JSS (in preparation).

EOF
}

cmd_validate() {
  echo "[validate] Looking for SAS at ${SAS_HOME} ..."
  if ! resolve_sas; then
    echo "[validate] ERROR: no SAS executable found under ${SAS_HOME}" >&2
    echo "[validate]        Mount your SAS install with -v /your/sas:/opt/sas:ro" >&2
    exit 2
  fi
  echo "[validate] Found SAS: ${SAS_BIN}"

  echo "[validate] Macros directory: ${GBTM_MACROS}"
  ls "${GBTM_MACROS}" || true

  # Tiny SAS smoke test
  TMP_SAS=$(mktemp --suffix=.sas)
  TMP_LOG=$(mktemp --suffix=.log)
  cat >"${TMP_SAS}" <<'SAS'
%put NOTE: Traj2 entrypoint validate -- SAS is reachable.;
proc options option=fullstimer; run;
SAS
  echo "[validate] Running SAS smoke test ..."
  if "${SAS_BIN}" -nodms -log "${TMP_LOG}" -sysin "${TMP_SAS}"; then
    echo "[validate] OK. SAS ran. Log:"
    tail -n 20 "${TMP_LOG}"
    rm -f "${TMP_SAS}" "${TMP_LOG}"
    exit 0
  else
    echo "[validate] ERROR: SAS returned non-zero. Log tail:" >&2
    tail -n 40 "${TMP_LOG}" >&2
    rm -f "${TMP_SAS}" "${TMP_LOG}"
    exit 3
  fi
}

cmd_list() {
  echo "Available macro source files in ${GBTM_MACROS}:"
  find "${GBTM_MACROS}" -maxdepth 2 -type f \
    \( -iname "*.sas" -o -iname "*.docx" -o -iname "*.rtf" \) \
    -printf "  %P\n" | sort
}

cmd_run() {
  local model="" nclass="" order="2" T="12"
  while [ $# -gt 0 ]; do
    case "$1" in
      --model=*)  model="${1#*=}"; shift ;;
      --model)    model="$2"; shift 2 ;;
      --nclass=*) nclass="${1#*=}"; shift ;;
      --nclass)   nclass="$2"; shift 2 ;;
      --order=*)  order="${1#*=}"; shift ;;
      --order)    order="$2"; shift 2 ;;
      --T=*)      T="${1#*=}"; shift ;;
      --T)        T="$2"; shift 2 ;;
      *) echo "Unknown arg: $1" >&2; exit 64 ;;
    esac
  done

  if [ -z "${model}" ] || [ -z "${nclass}" ]; then
    echo "Usage: run --model={ordinal|continuous|zip|zinb} --nclass=N [--order=2] [--T=12]" >&2
    exit 64
  fi

  if ! resolve_sas; then
    echo "ERROR: SAS not found at ${SAS_HOME}. Mount it with -v /your/sas:/opt/sas:ro" >&2
    exit 2
  fi

  case "${model}" in
    ordinal|continuous|zip|zinb) ;;
    *) echo "ERROR: --model must be one of: ordinal, continuous, zip, zinb" >&2; exit 64 ;;
  esac

  echo "[run] model=${model} nclass=${nclass} order=${order} T=${T}"
  echo "[run] data=${GBTM_DATA} output=${GBTM_OUTPUT}"
  echo "[run] SAS=${SAS_BIN}"

  # Generate a tiny driver SAS program. Real users will more likely supply
  # their own driver; this is a thin wrapper for the simplest case.
  local driver
  driver=$(mktemp --suffix=.sas)
  cat >"${driver}" <<SAS
options nofmterr;
libname GBTMOUT "${GBTM_OUTPUT}";
libname GBTMDAT "${GBTM_DATA}";
%let GBTM_MACROS = ${GBTM_MACROS};

/* User must place a driver named driver_<model>.sas under their data mount,
   or rely on the bundled macros_sas/run_<model>.sas if present. */
%macro _traj2_run;
  %let drv_user   = ${GBTM_DATA}/driver_${model}.sas;
  %let drv_bundle = ${GBTM_MACROS}/sas/run_${model}.sas;
  %if %sysfunc(fileexist(&drv_user.)) %then %do;
    %put NOTE: Using user driver &drv_user.;
    %include "&drv_user.";
  %end;
  %else %if %sysfunc(fileexist(&drv_bundle.)) %then %do;
    %put NOTE: Using bundled driver &drv_bundle.;
    %include "&drv_bundle.";
  %end;
  %else %do;
    %put ERROR: No driver found. Provide ${GBTM_DATA}/driver_${model}.sas
                or ship ${GBTM_MACROS}/sas/run_${model}.sas in the image.;
  %end;
%mend;
%_traj2_run;
SAS

  local logf="${GBTM_OUTPUT}/traj2_${model}_K${nclass}.log"
  local lstf="${GBTM_OUTPUT}/traj2_${model}_K${nclass}.lst"
  echo "[run] Log: ${logf}"
  echo "[run] LST: ${lstf}"

  "${SAS_BIN}" -nodms -log "${logf}" -print "${lstf}" -sysin "${driver}" \
    -set NCLASS "${nclass}" -set GBTM_ORDER "${order}" -set GBTM_T "${T}"
  rc=$?
  rm -f "${driver}"
  exit ${rc}
}

cmd_shell() {
  exec /bin/bash
}

# Dispatch
case "${1:-help}" in
  help|-h|--help)   cmd_help ;;
  validate)         cmd_validate ;;
  list)             cmd_list ;;
  run)              shift; cmd_run "$@" ;;
  shell|bash)       cmd_shell ;;
  *) echo "Unknown command: $1"; cmd_help; exit 64 ;;
esac

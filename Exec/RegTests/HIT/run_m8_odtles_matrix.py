#!/usr/bin/env python3
import argparse
import csv
import json
import os
import re
import shlex
import subprocess
import time
from pathlib import Path


EVENT_RATES = [0, 50, 100, 250, 500]
SUBSEGMENTS = [3, 5, 7]

RE_LES_STEPPED = re.compile(
    r"ODTLES: stepped (?P<stepped_entries>\d+) local lines over dt_LES with neutral SGS return\. "
    r"\[(?P<fields>.+)\]$"
)
RE_LES_REJECTED = re.compile(
    r"ODTLES: rejected (?P<rejected_total>\d+) owner-direction entries .* "
    r"\[(?P<fields>.+)\]$"
)
RE_FINAL_STEP = re.compile(
    r"STEP = (?P<step>\d+) TIME = (?P<time>[0-9eE+\-\.]+) DT = (?P<dt>[0-9eE+\-\.]+)"
)
RE_RUN_TIME = re.compile(r"Run time = (?P<runtime>[0-9eE+\-\.]+)")
RE_RUN_TIME_WO_INIT = re.compile(r"Run time w/o init = (?P<runtime>[0-9eE+\-\.]+)")
RE_TIME_REAL = re.compile(r"real (?P<real>[0-9eE+\-\.]+)")
RE_TIME_USER = re.compile(r"user (?P<user>[0-9eE+\-\.]+)")
RE_TIME_SYS = re.compile(r"sys (?P<sys>[0-9eE+\-\.]+)")


def parse_fields(field_blob: str):
    out = {}
    for chunk in field_blob.split(","):
        part = chunk.strip()
        if "=" not in part:
            continue
        k, v = part.split("=", 1)
        out[k.strip()] = v.strip()
    return out


def to_int(value: str):
    return int(float(value))


def to_float(value: str):
    return float(value)


def parse_log(log_path: Path):
    data = {
        "stepped_lines": 0,
        "stepped_entries_total": 0,
        "attempted_events_total": 0,
        "applied_events_total": 0,
        "rejected_state_model_events_total": 0,
        "rejected_events_total": 0,
        "diffusion_only_catchup_entries_total": 0,
        "directional_momentum_columns_total": 0,
        "rejected_support_total": 0,
        "rejected_boundary_total": 0,
        "rejected_amr_total": 0,
        "rejected_invalid_total": 0,
        "rejected_mixed_total": 0,
        "l1_umx_last": 0.0,
        "l1_umy_last": 0.0,
        "l1_umz_last": 0.0,
        "l1_ueden_last": 0.0,
        "l1_ueden_max_abs": 0.0,
        "last_step": None,
        "last_time": None,
        "last_dt": None,
        "run_time_reported": None,
        "run_time_wo_init_reported": None,
        "time_real": None,
        "time_user": None,
        "time_sys": None,
    }

    with log_path.open("r", errors="replace") as f:
        for line in f:
            s = line.strip()

            m = RE_LES_STEPPED.search(s)
            if m:
                fields = parse_fields(m.group("fields"))
                data["stepped_lines"] += 1
                data["stepped_entries_total"] += to_int(m.group("stepped_entries"))
                data["attempted_events_total"] += to_int(fields["attempted_events"])
                data["applied_events_total"] += to_int(fields["applied_events"])
                data["rejected_state_model_events_total"] += to_int(
                    fields["rejected_state_model_events"]
                )
                data["rejected_events_total"] += to_int(fields["rejected_events"])
                data["diffusion_only_catchup_entries_total"] += to_int(
                    fields["diffusion_only_catchup_entries"]
                )
                data["directional_momentum_columns_total"] += to_int(
                    fields["directional_momentum_columns"]
                )
                data["l1_umx_last"] = to_float(fields["l1_umx"])
                data["l1_umy_last"] = to_float(fields["l1_umy"])
                data["l1_umz_last"] = to_float(fields["l1_umz"])
                data["l1_ueden_last"] = to_float(fields["l1_ueden"])
                data["l1_ueden_max_abs"] = max(
                    data["l1_ueden_max_abs"], abs(data["l1_ueden_last"])
                )
                continue

            m = RE_LES_REJECTED.search(s)
            if m:
                fields = parse_fields(m.group("fields"))
                data["rejected_support_total"] += to_int(m.group("rejected_total"))
                data["rejected_boundary_total"] += to_int(fields["boundary_ghost"])
                data["rejected_amr_total"] += to_int(fields["amr_coarse_fine"])
                data["rejected_invalid_total"] += to_int(fields["invalid"])
                data["rejected_mixed_total"] += to_int(fields["mixed"])
                continue

            m = RE_FINAL_STEP.search(s)
            if m:
                data["last_step"] = int(m.group("step"))
                data["last_time"] = to_float(m.group("time"))
                data["last_dt"] = to_float(m.group("dt"))
                continue

            m = RE_RUN_TIME.search(s)
            if m:
                data["run_time_reported"] = to_float(m.group("runtime"))
                continue

            m = RE_RUN_TIME_WO_INIT.search(s)
            if m:
                data["run_time_wo_init_reported"] = to_float(m.group("runtime"))
                continue

            m = RE_TIME_REAL.search(s)
            if m:
                data["time_real"] = to_float(m.group("real"))
                continue

            m = RE_TIME_USER.search(s)
            if m:
                data["time_user"] = to_float(m.group("user"))
                continue

            m = RE_TIME_SYS.search(s)
            if m:
                data["time_sys"] = to_float(m.group("sys"))
                continue

    return data


def run_case(
    exe: Path,
    input_file: Path,
    ic_file: Path,
    out_dir: Path,
    event_rate: int,
    subsegments: int,
):
    case_name = f"er{event_rate:03d}_seg{subsegments}"
    case_dir = out_dir / case_name
    case_dir.mkdir(parents=True, exist_ok=True)
    case_input = case_dir / "input.inp"
    case_log = case_dir / "run.log"
    case_input.write_text(input_file.read_text())

    cmd = [
        "mpiexec",
        "-n",
        "32",
        str(exe),
        str(case_input.name),
        f"pelec.odt_event_rate={event_rate}",
        f"pelec.odt_subsegments_per_host_cell={subsegments}",
        f"prob.iname={ic_file}",
        "amr.plot_int=100000000",
        "amr.check_int=100000000",
        "amr.plot_file=plt",
        "amr.check_file=chk",
        "amr.data_log=datlog",
    ]

    with case_log.open("w") as logf:
        logf.write("# COMMAND\n")
        logf.write(" ".join(shlex.quote(x) for x in cmd) + "\n\n")
        logf.flush()
        t0 = time.perf_counter()
        proc = subprocess.run(
            ["/usr/bin/time", "-p"] + cmd,
            cwd=case_dir,
            stdout=logf,
            stderr=logf,
            check=False,
        )
        wall = time.perf_counter() - t0

    parsed = parse_log(case_log)
    parsed["event_rate"] = event_rate
    parsed["subsegments_per_host_cell"] = subsegments
    parsed["case_name"] = case_name
    parsed["exit_code"] = proc.returncode
    parsed["wall_time_seconds"] = wall
    parsed["command"] = " ".join(shlex.quote(x) for x in cmd)
    parsed["run_dir"] = str(case_dir.resolve())
    parsed["log_file"] = str(case_log.resolve())
    parsed["status"] = "ok" if proc.returncode == 0 else "failed"
    return parsed


def main():
    parser = argparse.ArgumentParser(
        description="Run M8 ODTLES HIT matrix on 32 MPI ranks and summarize diagnostics."
    )
    parser.add_argument(
        "--exe",
        required=True,
        help="Path to PeleC HIT executable (e.g. ./PeleC3d.gnu.MPI.ex).",
    )
    parser.add_argument(
        "--input",
        default="hit-les.inp",
        help="Baseline input file.",
    )
    parser.add_argument(
        "--ic-file",
        default="hit_ic_4_32.dat",
        help="Path to HIT initial condition file.",
    )
    parser.add_argument(
        "--outdir",
        default="m8_odtles_runs",
        help="Output directory for all cases and summaries.",
    )
    args = parser.parse_args()

    exe = Path(args.exe).resolve()
    input_file = Path(args.input).resolve()
    ic_file = Path(args.ic_file).resolve()
    out_dir = Path(args.outdir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    if not exe.exists():
        raise FileNotFoundError(f"Executable not found: {exe}")
    if not input_file.exists():
        raise FileNotFoundError(f"Input not found: {input_file}")
    if not ic_file.exists():
        raise FileNotFoundError(f"IC file not found: {ic_file}")

    summary = []
    for event_rate in EVENT_RATES:
        for subsegments in SUBSEGMENTS:
            print(
                f"[RUN] event_rate={event_rate:>3} subsegments={subsegments} (mpiexec -n 32)"
            )
            case = run_case(exe, input_file, ic_file, out_dir, event_rate, subsegments)
            print(
                f"      status={case['status']} exit={case['exit_code']} "
                f"wall={case['wall_time_seconds']:.3f}s attempted={case['attempted_events_total']} "
                f"applied={case['applied_events_total']} rejected={case['rejected_events_total']} "
                f"l1_ueden_max_abs={case['l1_ueden_max_abs']}"
            )
            summary.append(case)

    csv_fields = [
        "case_name",
        "event_rate",
        "subsegments_per_host_cell",
        "status",
        "exit_code",
        "wall_time_seconds",
        "run_time_reported",
        "run_time_wo_init_reported",
        "time_real",
        "time_user",
        "time_sys",
        "last_step",
        "last_time",
        "last_dt",
        "stepped_lines",
        "stepped_entries_total",
        "attempted_events_total",
        "applied_events_total",
        "rejected_state_model_events_total",
        "rejected_events_total",
        "diffusion_only_catchup_entries_total",
        "directional_momentum_columns_total",
        "rejected_support_total",
        "rejected_boundary_total",
        "rejected_amr_total",
        "rejected_invalid_total",
        "rejected_mixed_total",
        "l1_umx_last",
        "l1_umy_last",
        "l1_umz_last",
        "l1_ueden_last",
        "l1_ueden_max_abs",
        "run_dir",
        "log_file",
        "command",
    ]

    csv_path = out_dir / "m8_hit_odtles_summary.csv"
    with csv_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=csv_fields)
        writer.writeheader()
        for row in summary:
            writer.writerow(row)

    json_path = out_dir / "m8_hit_odtles_summary.json"
    with json_path.open("w") as f:
        json.dump(summary, f, indent=2)

    failed = [x for x in summary if x["status"] != "ok"]
    print(f"\nWrote: {csv_path}")
    print(f"Wrote: {json_path}")
    print(f"Completed {len(summary)} runs. Failed: {len(failed)}")
    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()

"""Sum MOOSE nonlinear iterations from a console log.

For each time step, MOOSE prints

    0 Nonlinear |R| = ...
    1 Nonlinear |R| = ...
    N Nonlinear |R| = ...
    Solve Converged!

Index 0 is the initial residual. The last index N is the nonlinear
iteration count for that step (the SNES iteration count at convergence).
Time Step 0 is the initial condition and has no solve.
"""

from __future__ import annotations

import re
from pathlib import Path

ANSI = re.compile(r"\x1b\[[0-9;]*m")
TIME_STEP = re.compile(r"^Time Step (\d+),")
NL = re.compile(r"^\s*(\d+) Nonlinear \|R\|")
SOLVER = re.compile(r"^  Solver Mode:\s+(.*\S)\s*$")


def parse(path: Path) -> dict:
    text = ANSI.sub("", path.read_text(errors="replace"))
    solver = None
    steps: list[dict] = []
    current: dict | None = None

    for line in text.splitlines():
        if solver is None:
            matched = SOLVER.match(line)
            if matched:
                solver = matched.group(1)
        step_match = TIME_STEP.match(line)
        if step_match:
            if current is not None:
                steps.append(current)
            current = {
                "step": int(step_match.group(1)),
                "nl": [],
                "converged": False,
                "failed": False,
            }
            continue
        if current is None:
            continue
        nl_match = NL.match(line)
        if nl_match:
            current["nl"].append(int(nl_match.group(1)))
            continue
        if "Solve Converged!" in line:
            current["converged"] = True
        elif "Solve Did Not Converge" in line or "Solve failed" in line.lower():
            current["failed"] = True

    if current is not None:
        steps.append(current)

    solved = [s for s in steps if s["nl"]]
    counts = []
    for s in solved:
        last = s["nl"][-1]
        if s["nl"] != list(range(last + 1)):
            raise SystemExit(f"{path.name}: step {s['step']} has non-sequential NL indices {s['nl']}")
        counts.append((s["step"], last, s["converged"], s["failed"]))

    return {
        "solver": solver,
        "n_solved_steps": len(counts),
        "n_converged": sum(1 for _, _, ok, _ in counts if ok),
        "n_failed": sum(1 for *_, failed in counts if failed),
        "nl_sum": sum(n for _, n, _, _ in counts),
        "nl_min": min(n for _, n, _, _ in counts),
        "nl_max": max(n for _, n, _, _ in counts),
        "counts": counts,
    }


def main() -> None:
    here = Path(__file__).resolve().parent
    logs = {
        "PJFNK": here / "s1_mohanty.log",
        "NEWTON": here / "step1_1.log",
    }
    results = {name: parse(path) for name, path in logs.items()}

    lines = ["case,time_step,nonlinear_iterations,converged"]
    for name, result in results.items():
        for step, n, converged, _failed in result["counts"]:
            lines.append(f"{name},{step},{n},{int(converged)}")
    (here / "nl_iterations_by_step.csv").write_text("\n".join(lines) + "\n")

    summary = ["case,solver_mode,solved_time_steps,converged,failed,nonlinear_iteration_sum,min_per_step,max_per_step"]
    for name, result in results.items():
        summary.append(
            f"{name},{result['solver']},{result['n_solved_steps']},{result['n_converged']},"
            f"{result['n_failed']},{result['nl_sum']},{result['nl_min']},{result['nl_max']}"
        )
    (here / "nl_iteration_sums.csv").write_text("\n".join(summary) + "\n")

    for name, result in results.items():
        print(
            f"{name}: solver={result['solver']}, steps={result['n_solved_steps']}, "
            f"converged={result['n_converged']}, failed={result['n_failed']}, "
            f"sum={result['nl_sum']}, per-step range={result['nl_min']}..{result['nl_max']}"
        )


if __name__ == "__main__":
    main()

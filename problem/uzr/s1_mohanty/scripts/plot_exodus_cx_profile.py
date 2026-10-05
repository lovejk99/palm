"""
从 MOOSE Exodus (.e) 提取多时刻位置–成分剖面 c(x)，右轴叠加温度 T(x)。

用法（须指定输入 .e；PNG 写在该 .e 所在目录）:
    python plot_exodus_cx_profile.py s1_mohanty_out.e
    python plot_exodus_cx_profile.py /path/to/case_out.e
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.io import netcdf_file

# =============================================================================
# 用户可调参数
# =============================================================================

# --- 时间采样 ---
# TARGET_TIMES 的数值单位，由 TIME_UNIT 指定（见下）
TARGET_TIMES: list[float] = [0.0, 5.0, 10.0, 20.0, 30.0]

# TARGET_TIMES 所用单位：'day'（天）或 'second'（秒，与 .e 中 time_whole 一致）
TIME_UNIT: str = "day"

# 当 TIME_UNIT='day' 时，将“天”换算为“秒”的因子，用于与 Exodus 的 time_whole 对齐。
# MOOSE 的 time_whole 始终以秒存储；你通常用“天”写 TARGET_TIMES，故需要 ×86400。
# 仅当 TIME_UNIT='second' 时本参数不参与换算（保留在此便于核对或特殊单位制）。
SECONDS_PER_DAY: float = 86400.0

# 时刻匹配：'nearest' 取最近帧；'exact' 要求偏差不超过 TIME_TOL
TIME_MATCH: str = "nearest"
TIME_TOL: float = 0.05  # 与 TIME_UNIT 同单位

# --- 空间范围（µm）---
X_MIN: float = 0.0
X_MAX: float = 300.0
X_LABEL: str = "Position (µm)"

# --- 场变量名（须与 Exodus name_nod_var 一致）---
VAR_CONCENTRATION: str = "c"
VAR_TEMPERATURE: str = "T"
C_AXIS_LABEL: str = "Zr mole fraction"
T_AXIS_LABEL: str = "Temperature (K)"

# --- 初始成分参考线 ---
SHOW_REF_LINE: bool = True
C0_INITIAL: float = 0.39

# 纵轴范围：None = 按数据自动留边距
C_YMIN: float | None = None
C_YMAX: float | None = None
T_YMIN: float | None = None
T_YMAX: float | None = None

# --- 温度曲线（当前算例 T 不随时间变，只画一条）---
SHOW_TEMPERATURE: bool = True
TEMPERATURE_TIME: float | None = None  # None = TARGET_TIMES 最后一项；单位同 TIME_UNIT

# --- 输出文件名（与输入 .e 同目录）---
# None = 由输入 .e 文件名推导（去掉 .e，并去掉末尾 _out，如 s1_mohanty_out.e → s1_mohanty.png）
OUTPUT_BASENAME: str | None = None
SAVE_PNG: bool = True

# --- 作图样式（matplotlib 默认）---
FIG_WIDTH_INCH: float = 8.0
FIG_HEIGHT_INCH: float = 5.0
SAVE_DPI: int = 200
LINEWIDTH: float = 1.5
SHOW_TITLE: bool = True
SHOW_ACTUAL_TIME_IN_LEGEND: bool = False  # True 时在图例显示匹配到的实际时间
LEGEND_LOC: str = "best"

# =============================================================================
# 实现
# =============================================================================


def cwd() -> Path:
    return Path(os.getcwd()).resolve()


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="从 Exodus (.e) 绘制多时刻 c(x) 剖面（右轴 T(x)）",
    )
    parser.add_argument(
        "exodus",
        type=str,
        help="输入 Exodus 文件路径（相对当前目录或绝对路径），例如 s1_mohanty_out.e",
    )
    return parser.parse_args(argv)


def resolve_exodus(exodus_arg: str) -> Path:
    path = Path(exodus_arg).expanduser()
    if not path.is_absolute():
        path = cwd() / path
    path = path.resolve()
    if not path.is_file():
        raise FileNotFoundError(f"未找到 Exodus 文件: {path}")
    return path


def output_basename(exodus_path: Path) -> str:
    if OUTPUT_BASENAME:
        return OUTPUT_BASENAME
    stem = exodus_path.stem
    if stem.endswith("_out"):
        stem = stem[: -len("_out")]
    return stem


def target_to_seconds(value: float) -> float:
    unit = TIME_UNIT.lower()
    if unit == "day":
        return value * SECONDS_PER_DAY
    if unit == "second":
        return value
    raise ValueError(f"TIME_UNIT 须为 'day' 或 'second'，当前为 {TIME_UNIT!r}")


def seconds_to_display(value_sec: float) -> tuple[float, str]:
    if TIME_UNIT.lower() == "day":
        return value_sec / SECONDS_PER_DAY, "d"
    return value_sec, "s"


def read_exodus(path: Path):
    with netcdf_file(path, mmap=False) as f:
        t = np.asarray(f.variables["time_whole"][:], dtype=float)
        x = np.asarray(f.variables["coordx"][:], dtype=float)
        names = [b"".join(n).decode().strip("\x00") for n in f.variables["name_nod_var"][:]]
        if VAR_CONCENTRATION not in names:
            raise KeyError(f"变量 '{VAR_CONCENTRATION}' 不在 {names}")
        if SHOW_TEMPERATURE and VAR_TEMPERATURE not in names:
            raise KeyError(f"变量 '{VAR_TEMPERATURE}' 不在 {names}")

        ic = names.index(VAR_CONCENTRATION)
        c_all = np.asarray(f.variables[f"vals_nod_var{ic + 1}"][:], dtype=float)
        T_all = None
        if SHOW_TEMPERATURE:
            iT = names.index(VAR_TEMPERATURE)
            T_all = np.asarray(f.variables[f"vals_nod_var{iT + 1}"][:], dtype=float)

    return t, x, c_all, T_all


def pick_time_index(t_sec: np.ndarray, target: float) -> tuple[int, float]:
    target_sec = target_to_seconds(target)
    idx = int(np.argmin(np.abs(t_sec - target_sec)))
    actual_sec = float(t_sec[idx])
    actual_disp, unit = seconds_to_display(actual_sec)
    target_disp, _ = seconds_to_display(target_sec)
    if TIME_MATCH == "exact" and abs(actual_disp - target_disp) > TIME_TOL:
        raise ValueError(
            f"未找到足够接近 {target_disp:g} {unit} 的时间步 "
            f"(最近 {actual_disp:.4f} {unit}, 容差 {TIME_TOL} {unit})"
        )
    return idx, actual_disp


def mask_position(x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    m = (x >= X_MIN) & (x <= X_MAX)
    if not np.any(m):
        raise ValueError(f"在 [{X_MIN}, {X_MAX}] µm 内无节点")
    order = np.argsort(x[m])
    return x[m][order], y[m][order]


def auto_ylim(data: np.ndarray, ymin, ymax, pad_frac: float = 0.06):
    if ymin is not None and ymax is not None:
        return ymin, ymax
    lo, hi = float(np.min(data)), float(np.max(data))
    span = hi - lo if hi > lo else max(abs(hi), 1.0) * 0.1
    if ymin is None and ymax is None:
        return lo - pad_frac * span, hi + pad_frac * span
    if ymin is None:
        return lo - pad_frac * span, ymax
    return ymin, hi + pad_frac * span


def plot_profiles(
    x: np.ndarray,
    profiles: list[tuple[float, float, np.ndarray]],
    T_x: np.ndarray | None,
    T_profile: np.ndarray | None,
    out_base: Path,
) -> list[Path]:
    fig, ax_c = plt.subplots(figsize=(FIG_WIDTH_INCH, FIG_HEIGHT_INCH))

    c_stack = []
    for target_t, actual_t, c_prof in profiles:
        if SHOW_ACTUAL_TIME_IN_LEGEND:
            label = f"{target_t:g} {TIME_UNIT} (actual {actual_t:.3g})"
        else:
            label = f"{target_t:g} d" if TIME_UNIT == "day" else f"{target_t:g} s"
        ax_c.plot(x, c_prof, lw=LINEWIDTH, label=label)
        c_stack.append(c_prof)

    if SHOW_REF_LINE:
        ax_c.axhline(
            C0_INITIAL,
            color="0.5",
            ls="--",
            lw=1.0,
            label=f"Initial ({C0_INITIAL:g})",
        )

    ax_c.set_xlabel(X_LABEL)
    ax_c.set_ylabel(C_AXIS_LABEL)
    c_ymin, c_ymax = auto_ylim(np.concatenate(c_stack), C_YMIN, C_YMAX)
    ax_c.set_ylim(c_ymin, c_ymax)
    ax_c.grid(True, alpha=0.3)

    if SHOW_TEMPERATURE and T_x is not None and T_profile is not None:
        ax_t = ax_c.twinx()
        ax_t.plot(
            T_x,
            T_profile,
            color="0.3",
            ls="--",
            lw=LINEWIDTH,
            label="T(x)",
        )
        ax_t.set_ylabel(T_AXIS_LABEL)
        t_ymin, t_ymax = auto_ylim(T_profile, T_YMIN, T_YMAX)
        ax_t.set_ylim(t_ymin, t_ymax)

        lines_c, labels_c = ax_c.get_legend_handles_labels()
        lines_t, labels_t = ax_t.get_legend_handles_labels()
        ax_c.legend(lines_c + lines_t, labels_c + labels_t, loc=LEGEND_LOC)
    else:
        ax_c.legend(loc=LEGEND_LOC)

    if SHOW_TITLE:
        ax_c.set_title("Zr concentration vs position")

    fig.tight_layout()

    saved: list[Path] = []
    if SAVE_PNG:
        p = out_base.with_suffix(".png")
        fig.savefig(p, dpi=SAVE_DPI, bbox_inches="tight")
        saved.append(p)
    plt.close(fig)
    return saved


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    print(f"[info] cwd: {cwd()}")

    exodus_path = resolve_exodus(args.exodus)
    print(f"[info] exodus: {exodus_path}")

    t_sec, x_raw, c_all, T_all = read_exodus(exodus_path)

    profiles: list[tuple[float, float, np.ndarray]] = []
    x_plot: np.ndarray | None = None
    for target_t in TARGET_TIMES:
        idx, actual_t = pick_time_index(t_sec, target_t)
        x_m, c_m = mask_position(x_raw, c_all[idx])
        if x_plot is None:
            x_plot = x_m
        profiles.append((target_t, actual_t, c_m))
        print(
            f"[info] target={target_t:g} {TIME_UNIT} -> "
            f"frame {idx}, actual={actual_t:.6f} {TIME_UNIT} ({t_sec[idx]:.0f} s)"
        )

    T_x, T_prof = None, None
    if SHOW_TEMPERATURE and T_all is not None:
        t_for_T = TEMPERATURE_TIME if TEMPERATURE_TIME is not None else TARGET_TIMES[-1]
        idx_t, actual_t = pick_time_index(t_sec, t_for_T)
        T_x, T_prof = mask_position(x_raw, T_all[idx_t])
        print(f"[info] T(x) from t={actual_t:.6f} {TIME_UNIT}")

    # 输出写在输入 .e 同目录，避免因启动目录不同而找不到图
    out_base = exodus_path.parent / output_basename(exodus_path)
    assert x_plot is not None
    saved = plot_profiles(x_plot, profiles, T_x, T_prof, out_base)

    if not saved:
        print("[warn] 未保存图片（SAVE_PNG = False）")
    else:
        for p in saved:
            print(f"[ok] {p}")


if __name__ == "__main__":
    try:
        main()
    except FileNotFoundError as exc:
        print(f"[error] {exc}", file=sys.stderr)
        sys.exit(1)

#!/usr/bin/env python3
"""sar-speckle-filtering — SAR斑点滤波

对 SAR（合成孔径雷达）强度影像执行自适应斑点噪声滤波。SAR 的相干成像
机制使影像叠加了乘性斑斑噪声（speckle），本 skill 实现三种经典方法：

- **Lee**（Lee 1980 最小均方误差滤波）：在滑窗内计算局部均值与局部方差，
  按 ``W = max(0, (var_local - var_noise) / var_local)`` 加权：
  ``out = mean + W * (pixel - mean)``。同质区域 W→0（输出局部均值，强平滑），
  边缘/纹理区 W→1（保留细节）。
- **Frost**（Frost 1982 指数衰减加权）：权重随到中心像元的距离按
  ``exp(-alpha * d)`` 衰减，衰减系数 ``alpha = 2 * var_noise / var_local``
  随局部方差自适应——方差小（同质）衰减快、平滑强，方差大（边缘）衰减慢、
  保细节。
- **multilook**（多视处理）：把影像切成 ``looks × looks`` 的块做平均，
  等效降低独立视数，是最朴素但稳健的降斑手段。

噪声模型假设为乘性：``I_obs = I_true * exp(N(0, sigma^2))``，因此局部
噪声方差取 ``var_noise = (noise_sigma * local_mean)^2``。

数据源：本地 SAR 强度/后向散射 GeoTIFF（线性幅度或强度），或使用
``--synthetic`` 生成一个平滑 σ⁰ 场叠加乘性斑斑噪声的模拟场景。

隐私声明 / Privacy：
- 默认完全离线运行，不发起任何网络请求。
- ``--synthetic`` 模式不读取任何外部数据。
- 所有处理在本地完成，不上传任何用户数据。

Usage:
    python sar-speckle-filtering.py --input sar.tif --filter lee --window 5
    python sar-speckle-filtering.py --bbox 116 39 117 40 --synthetic --filter frost

License: MIT
"""
from __future__ import annotations

import argparse
import datetime as _dt
import json
import os
import sys
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

VERSION = "1.0.0"
SKILL_NAME = "sar-speckle-filtering"

# ---- 复用共享核心库（本地 vendored，随脚本目录一起分发）----
try:
    from _geoskill_core.errors import (
        GeoSkillError, UsageError, ValidationError, ProcessError, to_exit_code,
    )
    from _geoskill_core.manifest import OutputManifest, OutputFile
except ImportError:  # pragma: no cover - fallback minimal definitions
    class GeoSkillError(Exception):
        def __init__(self, message: str, code: int = 7, kind: str = "EGeo", **kw):
            super().__init__(message)
            self.message, self.code, self.kind = message, code, kind

    class UsageError(GeoSkillError):
        def __init__(self, m, **k): super().__init__(m, code=2, kind="EUsage", **k)

    class ValidationError(GeoSkillError):
        def __init__(self, m, **k): super().__init__(m, code=6, kind="EValidate", **k)

    class ProcessError(GeoSkillError):
        def __init__(self, m, **k): super().__init__(m, code=7, kind="EProcess", **k)

    def to_exit_code(exc):
        return getattr(exc, "code", 7)

    OutputManifest = None
    OutputFile = None


FILTERS = ("lee", "frost", "multilook")


def _utc_now() -> str:
    return _dt.datetime.now(_dt.timezone.utc).isoformat().replace("+00:00", "Z")


# ---------------------------------------------------------------------------
# 核心算法
# ---------------------------------------------------------------------------
def lee_filter(
    img: np.ndarray,
    window: int = 5,
    noise_sigma: float = 0.3,
) -> np.ndarray:
    """Lee 自适应斑点滤波（Lee 1980 MMSE）。

    在 ``window × window`` 滑窗内计算局部均值/方差，按
    ``W = max(0, (var_local - var_noise)/var_local)`` 混合局部均值与原始像元。
    噪声方差用乘性模型 ``var_noise = (noise_sigma * local_mean)^2``。

    返回与输入同形的 float32 数组。
    """
    from scipy.ndimage import uniform_filter

    img = np.asarray(img, dtype=np.float32)
    if window < 3:
        raise UsageError(f"--window must be >= 3, got {window}", window=int(window))
    if window % 2 == 0:
        window += 1  # 保证奇数窗
    eps = 1e-12

    mean = uniform_filter(img, size=window, mode="reflect")
    mean_sq = uniform_filter(img * img, size=window, mode="reflect")
    var_local = np.maximum(mean_sq - mean * mean, 0.0)

    # 乘性噪声模型下的局部噪声方差
    var_noise = (float(noise_sigma) * mean) ** 2

    weight = np.zeros_like(var_local)
    denom = np.maximum(var_local, eps)
    mask = var_local > var_noise
    weight[mask] = (var_local[mask] - var_noise[mask]) / denom[mask]

    out = mean + weight * (img - mean)
    return out.astype(np.float32)


def frost_filter(
    img: np.ndarray,
    window: int = 5,
    noise_sigma: float = 0.3,
) -> np.ndarray:
    """Frost 指数衰减加权斑点滤波（Frost 1982）。

    局部衰减系数 ``alpha = 2 * var_noise / var_local``，权重
    ``w(d) = exp(-alpha * d)``（d 为到窗中心的欧氏距离），输出为窗内
    归一化加权平均。方差小→alpha 大→强平滑；方差大→alpha 小→保边缘。
    """
    from scipy.ndimage import uniform_filter

    img = np.asarray(img, dtype=np.float32)
    if window < 3:
        raise UsageError(f"--window must be >= 3, got {window}", window=int(window))
    if window % 2 == 0:
        window += 1
    half = window // 2
    eps = 1e-12

    mean = uniform_filter(img, size=window, mode="reflect")
    mean_sq = uniform_filter(img * img, size=window, mode="reflect")
    var_local = np.maximum(mean_sq - mean * mean, 0.0)
    var_noise = (float(noise_sigma) * mean) ** 2

    # alpha 随局部方差自适应；var_local 极小（同质）时给一个大衰减上限
    alpha = np.clip(2.0 * var_noise / np.maximum(var_local, eps), 0.0, 4.0)

    h, w = img.shape
    # 反射填充后按窗内偏移做向量化加权
    pad = np.pad(img, half, mode="reflect")
    num = np.zeros_like(img, dtype=np.float64)
    den = np.zeros_like(img, dtype=np.float64)
    for di in range(-half, half + 1):
        for dj in range(-half, half + 1):
            dist = float(np.hypot(di, dj))
            shifted = pad[half + di: half + di + h, half + dj: half + dj + w]
            weight = np.exp(-alpha * dist)
            num += weight * shifted
            den += weight
    out = num / np.maximum(den, eps)
    return out.astype(np.float32)


def multilook(img: np.ndarray, looks: int = 4) -> np.ndarray:
    """多视处理：把影像切成 ``looks × looks`` 块做平均。

    无法整除的边缘行列被裁掉。输出尺寸为 ``(H//looks, W//looks)``。
    """
    img = np.asarray(img, dtype=np.float32)
    if looks < 1:
        raise UsageError(f"--looks must be >= 1, got {looks}", looks=int(looks))
    if looks == 1:
        return img.copy()
    h, w = img.shape
    nh, nw = h // looks, w // looks
    if nh < 1 or nw < 1:
        raise ValidationError(
            f"image {h}x{w} too small for looks={looks}",
            height=int(h), width=int(w), looks=int(looks),
        )
    crop = img[: nh * looks, : nw * looks]
    blocked = crop.reshape(nh, looks, nw, looks)
    return blocked.mean(axis=(1, 3)).astype(np.float32)


def apply_filter(img: np.ndarray, method: str, window: int, looks: int,
                 noise_sigma: float) -> np.ndarray:
    """按方法分派滤波。"""
    if method == "lee":
        return lee_filter(img, window=window, noise_sigma=noise_sigma)
    if method == "frost":
        return frost_filter(img, window=window, noise_sigma=noise_sigma)
    if method == "multilook":
        return multilook(img, looks=looks)
    raise UsageError(
        f"unknown filter '{method}'. Choose from: {list(FILTERS)}",
        filter=method,
    )


# ---------------------------------------------------------------------------
# 合成数据：平滑 σ⁰ 场 + 乘性斑斑噪声
# ---------------------------------------------------------------------------
def generate_synthetic(
    bbox: List[float],
    width: int = 64,
    height: int = 64,
    noise_sigma: float = 0.3,
    seed: int = 42,
) -> Tuple[np.ndarray, Dict[str, Any]]:
    """生成一个 (1, H, W) 的含斑斑噪声 SAR 强度立方体。

    真值为一个平滑 σ⁰ 场（低频正弦基底 + 若干高斯地物），观测值在真值上
    叠加乘性斑斑噪声 ``exp(N(0, sigma^2))``。返回 (noisy_cube, info)，
    info 里带真值场用于质量验证。
    """
    rng = np.random.default_rng(seed)
    yy, xx = np.mgrid[0:height, 0:width].astype(np.float32)
    yn = yy / max(height - 1, 1)
    xn = xx / max(width - 1, 1)

    # 平滑真值 σ⁰ 场（线性强度，典型量级 0.02~0.2）
    truth = (
        0.10
        + 0.03 * np.sin(2.0 * np.pi * xn)
        + 0.03 * np.cos(2.0 * np.pi * yn)
    )
    # 叠加两个高斯"地物"（一亮一暗）
    truth += 0.08 * np.exp(-(((xn - 0.3) ** 2 + (yn - 0.3) ** 2) / 0.01))
    truth -= 0.05 * np.exp(-(((xn - 0.7) ** 2 + (yn - 0.7) ** 2) / 0.01))
    truth = np.clip(truth, 0.01, None).astype(np.float32)

    # 乘性斑斑噪声
    speckle = np.exp(rng.normal(0.0, noise_sigma, size=truth.shape)).astype(np.float32)
    noisy = (truth * speckle).astype(np.float32)

    info = {
        "bbox": bbox,
        "width": width,
        "height": height,
        "noise_sigma": float(noise_sigma),
        "truth_mean": float(np.mean(truth)),
        "truth_std": float(np.std(truth)),
        "noisy_std": float(np.std(noisy)),
        "truth_field": truth,
    }
    return noisy[np.newaxis, ...], info


# ---------------------------------------------------------------------------
# GeoTIFF I/O（直接用 rasterio，保证离线可用）
# ---------------------------------------------------------------------------
def write_geotiff(
    path: str,
    cube: np.ndarray,
    bbox: List[float],
    nodata: float = -9999.0,
) -> None:
    import rasterio
    from rasterio.transform import from_bounds

    if cube.ndim == 2:
        cube = cube[np.newaxis, ...]
    nb, h, w = cube.shape
    transform = from_bounds(bbox[0], bbox[1], bbox[2], bbox[3], w, h)
    profile = {
        "driver": "GTiff", "height": h, "width": w, "count": nb,
        "dtype": "float32", "crs": "EPSG:4326", "transform": transform,
        "nodata": nodata, "compress": "deflate",
    }
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with rasterio.open(path, "w", **profile) as dst:
        for b in range(nb):
            dst.write(cube[b].astype("float32"), b + 1)


def read_geotiff(path: str) -> Tuple[np.ndarray, List[float]]:
    import rasterio
    if not os.path.exists(path):
        raise UsageError(f"input raster not found: {path}", path=path)
    with rasterio.open(path) as src:
        cube = src.read().astype(np.float32)
        b = src.bounds
        bbox = [b.left, b.bottom, b.right, b.top]
    return cube, bbox


def read_geotiff_full(path: str) -> Tuple[np.ndarray, List[float], float]:
    """Read multi-band GeoTIFF, replace nodata with NaN, validate n_valid_pixels.

    Returns (cube_with_nan, bbox, nodata). Raises ValidationError if all pixels
    are NoData across all bands. nodata may be None if file has no nodata tag.
    """
    import rasterio
    if not os.path.exists(path):
        raise UsageError(f"input raster not found: {path}", path=path)
    with rasterio.open(path) as src:
        cube = src.read().astype(np.float32)
        b = src.bounds
        bbox = [b.left, b.bottom, b.right, b.top]
        nodata = src.nodata
    if nodata is not None and np.isfinite(nodata):
        cube = np.where(cube == nodata, np.nan, cube)
    n_valid = int(np.sum(np.isfinite(cube)))
    if n_valid == 0:
        raise ValidationError(
            f"input raster has no valid pixels (all {cube.size} are NoData={nodata})"
        )
    return cube, bbox, nodata


def validate_bbox(bbox: List[float]) -> None:
    """Validate bbox = [W, S, E, N]. Raise ValidationError on W>=E, S>=N, out-of-range,
    or cross-180° antipodal bbox."""
    if not (isinstance(bbox, (list, tuple)) and len(bbox) == 4):
        raise ValidationError(f"bbox must be 4 floats [W S E N], got {bbox}")
    w, s, e, n = [float(v) for v in bbox]
    if not (-180.0 <= w <= 180.0 and -180.0 <= e <= 180.0):
        raise ValidationError(
            f"longitude out of range [-180,180]: W={w}, E={e}"
        )
    if not (-90.0 <= s <= 90.0 and -90.0 <= n <= 90.0):
        raise ValidationError(
            f"latitude out of range [-90,90]: S={s}, N={n}"
        )
    if w >= e:
        if abs(e - (-180.0)) < 1e-9 and w > 0:
            raise ValidationError(
                f"cross-180° bbox not supported (W={w}, E={e}); "
                f"split into two non-antipodal bboxes"
            )
        raise ValidationError(f"W must be < E, got W={w}, E={e}")
    if s >= n:
        raise ValidationError(f"S must be < N, got S={s}, N={n}")
    if (e - w) < 0.001 or (n - s) < 0.001:
        raise ValidationError(
            f"bbox too small (<0.001°), got W={w},S={s},E={e},N={n}"
        )


# ---------------------------------------------------------------------------
# Manifest
# ---------------------------------------------------------------------------
def write_manifest(
    output_dir: str,
    args: argparse.Namespace,
    outputs: List[Dict[str, Any]],
    qa: Dict[str, Any],
    started_at: str,
    exit_code: int,
    bbox: List[float],
) -> Optional[str]:
    if OutputManifest is None:
        return None
    cmd = " ".join([SKILL_NAME] + sys.argv[1:])
    man = OutputManifest(
        skill=SKILL_NAME,
        skill_version=VERSION,
        command=cmd,
        started_at=started_at,
        finished_at=_utc_now(),
        exit_code=exit_code,
        inputs={
            "input": getattr(args, "input", None),
            "filter": getattr(args, "filter", None),
            "window": getattr(args, "window", None),
            "looks": getattr(args, "looks", None),
            "synthetic": bool(getattr(args, "synthetic", False)),
        },
        outputs=[OutputFile(**o) for o in outputs],
        qa=qa,
        software={"python": sys.version.split()[0], "skill": SKILL_NAME},
    )
    path = os.path.join(output_dir, "output-manifest.json")
    with open(path, "w", encoding="utf-8") as f:
        json.dump(man.to_dict(), f, ensure_ascii=False, indent=2, default=str)
    return path


# ---------------------------------------------------------------------------
# 主流程
# ---------------------------------------------------------------------------
def process(args: argparse.Namespace) -> int:
    started_at = _utc_now()
    output_dir = args.output_dir

    bbox = list(args.bbox) if args.bbox else None

    # 1) 获取数据
    synth_info: Optional[Dict[str, Any]] = None
    input_nodata: Optional[float] = None
    n_valid_pixels: Optional[int] = None
    if args.input and not args.synthetic:
        if bbox is not None:
            validate_bbox(bbox)
        cube, file_bbox, input_nodata = read_geotiff_full(args.input)
        bbox = bbox if bbox is not None else file_bbox
        n_valid_pixels = int(np.sum(np.isfinite(cube)))
        source_note = args.input
    else:
        if bbox is None:
            raise UsageError("provide --bbox (synthetic mode) or --input <raster>")
        validate_bbox(bbox)
        cube, synth_info = generate_synthetic(
            bbox, noise_sigma=args.noise_sigma,
        )
        n_valid_pixels = int(cube.size)
        source_note = "synthetic"

    if cube.size == 0:
        raise ValidationError("input raster is empty")
    if not np.all(np.isfinite(cube)):
        raise ValidationError("input raster contains non-finite values")

    # Now safe to create output dir
    os.makedirs(output_dir, exist_ok=True)

    # 2) 逐波段滤波
    filtered_bands = [
        apply_filter(cube[b], args.filter, args.window, args.looks, args.noise_sigma)
        for b in range(cube.shape[0])
    ]
    filtered = np.stack(filtered_bands, axis=0)

    # 3) 写出产物
    out_tif = os.path.join(output_dir, "filtered.tif")
    write_geotiff(out_tif, filtered, bbox)

    # QA：逐波段方差下降比 / 均值保持
    in_std = [float(np.std(cube[b])) for b in range(cube.shape[0])]
    out_std = [float(np.std(filtered_bands[b])) for b in range(len(filtered_bands))]
    in_mean = [float(np.mean(cube[b])) for b in range(cube.shape[0])]
    out_mean = [float(np.mean(filtered_bands[b])) for b in range(len(filtered_bands))]
    std_reduction = [
        (in_std[b] - out_std[b]) / in_std[b] if in_std[b] > 0 else 0.0
        for b in range(cube.shape[0])
    ]

    params = {
        "filter": args.filter,
        "window": int(args.window),
        "looks": int(args.looks),
        "noise_sigma": float(args.noise_sigma),
        "source": source_note,
        "input_shape": list(cube.shape),
        "output_shape": list(filtered.shape),
        "std_per_band_input": in_std,
        "std_per_band_output": out_std,
        "std_reduction_ratio": std_reduction,
        "mean_per_band_input": in_mean,
        "mean_per_band_output": out_mean,
    }
    params_path = os.path.join(output_dir, "filter_params.json")
    with open(params_path, "w", encoding="utf-8") as f:
        json.dump(params, f, ensure_ascii=False, indent=2)

    qa: Dict[str, Any] = {
        "source": source_note,
        "filter": args.filter,
        "n_bands": int(cube.shape[0]),
        "std_reduction_ratio_band0": std_reduction[0],
        "mean_preservation_band0": (
            out_mean[0] / in_mean[0] if in_mean[0] != 0 else 0.0
        ),
        "n_valid_pixels": n_valid_pixels,
        "input_nodata": input_nodata,
    }
    if synth_info is not None:
        qa["synthetic_truth_std"] = synth_info["truth_std"]
        qa["synthetic_noisy_std"] = synth_info["noisy_std"]

    outputs = [
        {"path": out_tif, "kind": "raster", "crs_epsg": 4326,
         "bbox_wgs84": bbox, "band_count": int(filtered.shape[0])},
        {"path": params_path, "kind": "json"},
    ]

    man_path = write_manifest(output_dir, args, outputs, qa, started_at, 0, bbox)

    if not args.quiet:
        print(f"[{SKILL_NAME}] source: {source_note}")
        print(f"[{SKILL_NAME}] filter: {args.filter}  window: {args.window}  looks: {args.looks}")
        print(f"[{SKILL_NAME}] input shape:  {cube.shape}")
        print(f"[{SKILL_NAME}] output shape: {filtered.shape}")
        print(f"[{SKILL_NAME}] output: {out_tif}")
        print(f"[{SKILL_NAME}] params: {params_path}")
        if man_path:
            print(f"[{SKILL_NAME}] manifest: {man_path}")
        print(f"[{SKILL_NAME}] band0 std reduction: {std_reduction[0]:.3f}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog=SKILL_NAME,
        description="SAR speckle filtering (Lee / Frost / multilook) for intensity imagery.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--bbox", nargs=4, type=float, metavar=("W", "S", "E", "N"),
                   help="geographic extent [minLon minLat maxLon maxLat]")
    p.add_argument("--input", help="input SAR intensity/backscatter GeoTIFF")
    p.add_argument("--filter", default="lee", choices=list(FILTERS),
                   help="speckle filter method (default: lee)")
    p.add_argument("--window", type=int, default=5,
                   help="filter window size in pixels for lee/frost (default: 5)")
    p.add_argument("--looks", type=int, default=4,
                   help="number of looks (block size) for multilook (default: 4)")
    p.add_argument("--noise-sigma", type=float, default=0.3,
                   help="relative multiplicative noise std (default: 0.3)")
    p.add_argument("--synthetic", action="store_true",
                   help="generate a physics-consistent synthetic scene (offline)")
    p.add_argument("--output-dir", default="./output", help="output directory")
    p.add_argument("--quiet", action="store_true", help="suppress console output")
    p.add_argument("--version", action="version", version=f"%(prog)s {VERSION}")
    return p


def main(argv: Optional[List[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return process(args)
    except GeoSkillError as exc:
        print(f"[{SKILL_NAME}] ERROR [{exc.kind}] {exc.message}", file=sys.stderr)
        return exc.code
    except KeyboardInterrupt:
        return 130
    except Exception as exc:  # noqa: BLE001
        print(f"[{SKILL_NAME}] ERROR {exc}", file=sys.stderr)
        return to_exit_code(exc)


if __name__ == "__main__":
    sys.exit(main())

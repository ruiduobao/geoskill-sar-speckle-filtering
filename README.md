# SAR Speckle Filtering (geoskill-sar-speckle-filtering)

> Applies Lee/Frost/multilook speckle filtering to SAR imagery to suppress coherent speckle noise while preserving edges. Lee/Frost/multilook speckle filtering for SAR imagery. Takes single-band or multi-band SAR intensity GeoTIFF input (or generates a scene with multiplicative speckle noise via --synthetic), and outputs the filtered GeoTIFF + parameter JSON.

---

## 1. Overview

(Fill in 2-3 paragraphs of introduction here: features, application scenarios, core algorithms.)

## 3. Quick Start

```bash
pip install -r requirements.txt
python geoskill-sar-speckle-filtering.py --bbox 116 39 117 40 --synthetic --output-dir ./out
```

## 4. CLI Parameters

Run `python <skill>.py --help` for the full list. Common parameters:

| Parameter | Type | Description |
|---|---|---|
| `--bbox` | `float[4]` | WGS84 bounding box `min_lon min_lat max_lon max_lat` |
| `--input` | `path` | Local input file (GeoJSON/GeoTIFF/etc.) |
| `--output-dir` | `path` | Output directory (default `./output`) |
| `--synthetic` | `flag` | Use synthetic data instead of real input |
| `--quiet` | `flag` | Suppress non-essential stdout |

## 5. Input / Output

| File | Format | Description |
|---|---|---|
| `result.tif` | GeoTIFF | Primary output |
| `output-manifest.json` | JSON | Run manifest |

## 6. Technical Principle

(see SKILL.md for details)

## 8. License

MIT License. See [`LICENSE`](./LICENSE) for full text.

---

> 中文完整版：[references/SKILL.zh-CN.md](references/SKILL.zh-CN.md)

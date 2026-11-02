---
name: geoskill-sar-speckle-filtering
description: '对 SAR 影像执行 Lee/Frost/多视 斑点噪声滤波，抑制相干斑噪声并保持边缘。Lee/Frost/multilook speckle filtering for SAR imagery. 输入单波段或多波段 SAR 强度 GeoTIFF（或用 --synthetic 生成含乘性斑斑噪声的场景），输出滤波后 GeoTIFF + 参数 JSON。'
---

# SAR斑点滤波 | SAR Speckle Filtering

(To be filled in: 2–3 paragraphs introducing the functionality, application scenarios, and core algorithm.)

## Dependencies / 依赖

```bash
pip install 'numpy' 'rasterio' 'scipy'
```

## Usage / 使用方法

### Basic Usage

```bash
python geoskill-sar-speckle-filtering.py --bbox 116.0 39.0 117.0 40.0 [other parameters]
```

### Example 1 (Synthetic Data, Offline)

```bash
python geoskill-sar-speckle-filtering.py --bbox 116.0 39.0 117.0 40.0 --synthetic --output-dir ./out
```

### Examples 2-5

(Add at least 4 real usage examples.)

## Output / 输出

| File | Format | Description |
|---|---|---|
| `result.tif` | GeoTIFF | Primary output |
| `output-manifest.json` | JSON | Run manifest |

## Data Source / 数据源 / Source

(Specify the data source: free satellite data / local input / synthetic.)

## Privacy / 隐私声明 / Privacy

- Runs offline by default; `--synthetic` mode requires no network at all.
- All processing is done locally; no user data is uploaded.

## License / License

MIT

---

> 中文完整版：需要时再读 [references/SKILL.zh-CN.md](references/SKILL.zh-CN.md)，不要默认全量加载。

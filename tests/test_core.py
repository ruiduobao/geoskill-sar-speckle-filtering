"""Core algorithm tests for sar-speckle-filtering."""
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(__file__))
from conftest import mod as ssf


def _speckled(seed=0, shape=(64, 64), sigma=0.3, base=0.1):
    """恒定真值场 + 乘性斑斑噪声。"""
    rng = np.random.default_rng(seed)
    truth = np.full(shape, base, dtype=np.float32)
    noisy = truth * np.exp(rng.normal(0.0, sigma, size=shape)).astype(np.float32)
    return truth, noisy.astype(np.float32)


class TestLee:
    def test_reduces_variance(self):
        truth, noisy = _speckled(sigma=0.3)
        out = ssf.lee_filter(noisy, window=5, noise_sigma=0.3)
        assert out.shape == noisy.shape
        # 滤波后方差应明显下降
        assert np.std(out) < np.std(noisy) * 0.8

    def test_preserves_mean(self):
        truth, noisy = _speckled(sigma=0.3)
        out = ssf.lee_filter(noisy, window=5, noise_sigma=0.3)
        # 均值基本保持（乘性噪声均值≈真值）
        assert abs(np.mean(out) - np.mean(truth)) / np.mean(truth) < 0.15

    def test_constant_image_unchanged(self):
        img = np.full((32, 32), 0.5, dtype=np.float32)
        out = ssf.lee_filter(img, window=5, noise_sigma=0.3)
        # 常数场局部方差=0 → W=0 → 输出=均值=常数
        np.testing.assert_allclose(out, img, atol=1e-5)

    def test_window_too_small_raises(self):
        img = np.ones((16, 16), dtype=np.float32)
        with pytest.raises(ssf.UsageError):
            ssf.lee_filter(img, window=1)


class TestFrost:
    def test_reduces_variance(self):
        truth, noisy = _speckled(sigma=0.3)
        out = ssf.frost_filter(noisy, window=5, noise_sigma=0.3)
        assert out.shape == noisy.shape
        assert np.std(out) < np.std(noisy)

    def test_preserves_mean(self):
        truth, noisy = _speckled(sigma=0.3)
        out = ssf.frost_filter(noisy, window=5, noise_sigma=0.3)
        assert abs(np.mean(out) - np.mean(truth)) / np.mean(truth) < 0.2

    def test_constant_image_unchanged(self):
        img = np.full((32, 32), 0.4, dtype=np.float32)
        out = ssf.frost_filter(img, window=5, noise_sigma=0.3)
        np.testing.assert_allclose(out, img, atol=1e-4)


class TestMultilook:
    def test_downsamples_by_looks(self):
        img = np.random.default_rng(1).uniform(0.05, 0.2, (64, 64)).astype(np.float32)
        out = ssf.multilook(img, looks=4)
        assert out.shape == (16, 16)

    def test_reduces_variance(self):
        truth, noisy = _speckled(sigma=0.3)
        out = ssf.multilook(noisy, looks=4)
        assert np.std(out) < np.std(noisy)

    def test_looks_one_identity(self):
        img = np.random.default_rng(2).uniform(0, 1, (8, 8)).astype(np.float32)
        out = ssf.multilook(img, looks=1)
        np.testing.assert_allclose(out, img)

    def test_too_small_raises(self):
        img = np.ones((3, 3), dtype=np.float32)
        with pytest.raises(ssf.ValidationError):
            ssf.multilook(img, looks=4)

    def test_block_average_exact(self):
        img = np.array([[1, 3], [5, 7]], dtype=np.float32)
        out = ssf.multilook(img, looks=2)
        assert out.shape == (1, 1)
        np.testing.assert_allclose(out[0, 0], 4.0)  # mean(1,3,5,7)=4


class TestApplyFilter:
    def test_dispatch_unknown_raises(self):
        img = np.ones((8, 8), dtype=np.float32)
        with pytest.raises(ssf.UsageError):
            ssf.apply_filter(img, "bilateral", 5, 4, 0.3)


class TestSynthetic:
    def test_shape_and_info(self):
        cube, info = ssf.generate_synthetic([116, 39, 117, 40])
        assert cube.ndim == 3
        assert cube.shape == (1, 64, 64)
        assert info["truth_std"] >= 0
        assert "truth_field" in info
        assert np.all(cube > 0)

    def test_noisy_more_variable_than_truth(self):
        cube, info = ssf.generate_synthetic([116, 39, 117, 40], noise_sigma=0.4)
        truth = info["truth_field"]
        # 斑斑噪声应显著增大影像方差
        assert np.std(cube[0]) > np.std(truth)

    def test_filter_recovers_toward_truth(self):
        """Lee 滤波应使结果比噪声影像更接近真值场。"""
        cube, info = ssf.generate_synthetic([116, 39, 117, 40], noise_sigma=0.35, seed=7)
        truth = info["truth_field"]
        noisy = cube[0]
        filtered = ssf.lee_filter(noisy, window=5, noise_sigma=0.35)
        err_noisy = np.mean((noisy - truth) ** 2)
        err_filt = np.mean((filtered - truth) ** 2)
        assert err_filt < err_noisy


class TestGeoTiffIO:
    def test_write_and_read_roundtrip(self, tmp_path):
        cube = np.random.default_rng(3).uniform(0.01, 0.3, (2, 16, 16)).astype(np.float32)
        bbox = [116.0, 39.0, 117.0, 40.0]
        path = str(tmp_path / "test.tif")
        ssf.write_geotiff(path, cube, bbox)
        assert os.path.exists(path)
        back, rbbox = ssf.read_geotiff(path)
        assert back.shape == cube.shape
        np.testing.assert_allclose(rbbox, bbox, atol=1e-6)
        np.testing.assert_allclose(back, cube, atol=1e-5)

    def test_read_missing_file_raises(self):
        with pytest.raises(ssf.UsageError):
            ssf.read_geotiff("/nonexistent/path/file.tif")

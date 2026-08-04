"""CLI argument parsing tests for sar-speckle-filtering."""
import subprocess
import sys
import os

sys.path.insert(0, os.path.dirname(__file__))
from conftest import mod, SCRIPT


def run_cli(args, timeout=120):
    return subprocess.run(
        [sys.executable, SCRIPT] + args,
        capture_output=True, text=True, timeout=timeout,
    )


class TestCLIBasics:
    def test_help(self):
        r = run_cli(["--help"])
        assert r.returncode == 0
        assert "--bbox" in r.stdout
        assert "--synthetic" in r.stdout
        assert "--filter" in r.stdout
        assert "--window" in r.stdout

    def test_version(self):
        r = run_cli(["--version"])
        assert r.returncode == 0
        assert "1.0.0" in r.stdout

    def test_no_args_exit_2(self):
        """No arguments → UsageError → exit 2."""
        r = run_cli([])
        assert r.returncode == 2

    def test_bad_filter_exit_2(self):
        """Invalid filter choice → argparse error → exit 2."""
        r = run_cli(["--bbox", "116", "39", "117", "40", "--filter", "bad"])
        assert r.returncode == 2

    def test_input_file_not_found_exit_2(self):
        """Non-existent input file → UsageError → exit 2."""
        r = run_cli(["--input", "nonexistent.tif"])
        assert r.returncode == 2

    def test_synthetic_without_bbox_exit_2(self):
        """--synthetic without --bbox → UsageError → exit 2."""
        r = run_cli(["--synthetic"])
        assert r.returncode == 2


class TestCLIRun:
    def test_synthetic_lee_runs(self, tmp_path):
        out = str(tmp_path / "out_lee")
        r = run_cli([
            "--bbox", "116", "39", "117", "40",
            "--synthetic", "--filter", "lee", "--window", "5",
            "--output-dir", out, "--quiet",
        ])
        assert r.returncode == 0, r.stderr
        assert os.path.exists(os.path.join(out, "filtered.tif"))
        assert os.path.exists(os.path.join(out, "output-manifest.json"))

    def test_synthetic_frost_runs(self, tmp_path):
        out = str(tmp_path / "out_frost")
        r = run_cli([
            "--bbox", "116", "39", "117", "40",
            "--synthetic", "--filter", "frost", "--quiet",
            "--output-dir", out,
        ])
        assert r.returncode == 0, r.stderr
        assert os.path.exists(os.path.join(out, "filtered.tif"))

    def test_synthetic_multilook_runs(self, tmp_path):
        out = str(tmp_path / "out_ml")
        r = run_cli([
            "--bbox", "116", "39", "117", "40",
            "--synthetic", "--filter", "multilook", "--looks", "4",
            "--output-dir", out, "--quiet",
        ])
        assert r.returncode == 0, r.stderr
        assert os.path.exists(os.path.join(out, "filtered.tif"))

    def test_bbox_only_auto_synthetic(self, tmp_path):
        """Only --bbox (no --input, no --synthetic) → auto synthetic → rc 0."""
        out = str(tmp_path / "out_auto")
        r = run_cli([
            "--bbox", "116.39", "39.90", "116.40", "39.91",
            "--output-dir", out, "--quiet",
        ])
        assert r.returncode == 0, r.stderr
        assert os.path.exists(os.path.join(out, "output-manifest.json"))

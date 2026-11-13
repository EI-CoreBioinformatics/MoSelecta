from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

from moselecta.compare import compare_file_pair


def package_root() -> Path:
    return Path(__file__).resolve().parents[1]


def run_command(args: list[str]) -> subprocess.CompletedProcess[str]:
    root = package_root()
    env = os.environ.copy()
    env["PYTHONPATH"] = str(root / "src")
    return subprocess.run(
        args,
        cwd=root,
        env=env,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=True,
    )


def test_compare_tsv_tolerates_last_decimal_numeric_difference(tmp_path: Path) -> None:
    left = tmp_path / "left.tsv"
    right = tmp_path / "right.tsv"
    left.write_text("id\tscore\nmodel1\t0.810938\n")
    right.write_text("id\tscore\nmodel1\t0.810937\n")

    result = compare_file_pair(left, right, abs_tol=1e-6, rel_tol=1e-9)
    assert result.status == "equivalent"
    assert result.equivalent


def test_compare_tsv_rejects_string_difference(tmp_path: Path) -> None:
    left = tmp_path / "left.tsv"
    right = tmp_path / "right.tsv"
    left.write_text("id\tscore\nmodel1\t0.5\n")
    right.write_text("id\tscore\nmodel2\t0.5\n")

    result = compare_file_pair(left, right, abs_tol=1e-6, rel_tol=1e-9)
    assert result.status == "different"
    assert not result.equivalent



def test_compare_tsv_tolerates_numeric_difference_inside_component_string(tmp_path: Path) -> None:
    left = tmp_path / "left.tsv"
    right = tmp_path / "right.tsv"
    left.write_text("id\tcomponents\nmodel1\tAnnevo:0.701416;Helixer:0.480938\n")
    right.write_text("id\tcomponents\nmodel1\tAnnevo:0.701416;Helixer:0.480937\n")

    result = compare_file_pair(left, right, abs_tol=1e-6, rel_tol=1e-9)
    assert result.status == "equivalent"
    assert result.equivalent


def test_compare_tsv_rejects_component_string_template_difference(tmp_path: Path) -> None:
    left = tmp_path / "left.tsv"
    right = tmp_path / "right.tsv"
    left.write_text("id\tcomponents\nmodel1\tAnnevo:0.701416;Helixer:0.480938\n")
    right.write_text("id\tcomponents\nmodel1\tAnnevo:0.701416;Tiberius:0.480937\n")

    result = compare_file_pair(left, right, abs_tol=1e-6, rel_tol=1e-9)
    assert result.status == "different"
    assert not result.equivalent


def test_compare_cli_smoke(tmp_path: Path) -> None:
    left = tmp_path / "left.tsv"
    right = tmp_path / "right.tsv"
    left.write_text("id\tscore\nmodel1\t0.810938\n")
    right.write_text("id\tscore\nmodel1\t0.810937\n")
    result = run_command([
        sys.executable,
        "-m",
        "moselecta.compare",
        "--files",
        str(left),
        str(right),
    ])
    assert result.returncode == 0
    assert "EQUIVALENT" in result.stdout

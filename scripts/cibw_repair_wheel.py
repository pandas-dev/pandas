#!/usr/bin/env python3
"""
Repair wheel and inject SBOM for PEP 770 compliance.

This script:
1. Runs the platform-specific wheel repair tool (auditwheel/delvewheel)
2. Injects the vendored code SBOM into the repaired wheel

Usage:
    python scripts/cibw_repair_wheel.py <wheel> <dest_dir>
"""

import argparse
import base64
import hashlib
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import zipfile


def get_wheel_dist_info(wheel_path: Path) -> str:
    """Get the .dist-info directory name from a wheel."""
    with zipfile.ZipFile(wheel_path, "r") as zf:
        for name in zf.namelist():
            if ".dist-info/" in name:
                # Extract just the dist-info directory name
                return name.split("/")[0]
    raise ValueError(f"No .dist-info directory found in {wheel_path}")


def inject_sbom(wheel_path: Path, sbom_path: Path) -> None:
    """Inject SBOM into wheel's .dist-info/sboms/ directory."""
    dist_info = get_wheel_dist_info(wheel_path)
    sbom_wheel_path = f"{dist_info}/sboms/{sbom_path.name}"

    # Read existing wheel contents
    with zipfile.ZipFile(wheel_path, "a") as zf:
        # Check if SBOM already exists (e.g., from auditwheel)
        existing_sboms = [n for n in zf.namelist() if "/sboms/" in n]
        if existing_sboms:
            print(f"  Existing SBOMs in wheel: {existing_sboms}")

        # Add our vendored code SBOM
        print(f"  Adding {sbom_wheel_path}")
        zf.write(sbom_path, sbom_wheel_path)

    # Update RECORD file
    update_record(wheel_path, sbom_wheel_path, sbom_path)


def update_record(wheel_path: Path, sbom_wheel_path: str, sbom_path: Path) -> None:
    """Update the RECORD file in the wheel to include the SBOM."""
    # Calculate hash of SBOM file
    with open(sbom_path, "rb") as f:
        content = f.read()
        sha256_hash = hashlib.sha256(content).digest()
        hash_digest = base64.urlsafe_b64encode(sha256_hash).rstrip(b"=").decode("ascii")

    # Format: path,hash,size
    record_line = f"{sbom_wheel_path},sha256={hash_digest},{len(content)}"

    # Read existing RECORD, append new entry
    with zipfile.ZipFile(wheel_path, "r") as zf:
        dist_info = get_wheel_dist_info(wheel_path)
        record_path = f"{dist_info}/RECORD"
        record_content = zf.read(record_path).decode("utf-8")

    # Append SBOM entry to RECORD
    record_lines = record_content.rstrip("\n").split("\n")
    record_lines.append(record_line)
    new_record = "\n".join(record_lines) + "\n"

    # Rewrite wheel with updated RECORD
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_path = Path(tmpdir)

        # Extract wheel
        with zipfile.ZipFile(wheel_path, "r") as zf:
            zf.extractall(tmp_path)

        # Update RECORD
        (tmp_path / record_path).write_text(new_record, encoding="utf-8")

        # Repack wheel
        wheel_path.unlink()
        with zipfile.ZipFile(wheel_path, "w", compression=zipfile.ZIP_DEFLATED) as zf:
            for file_path in tmp_path.rglob("*"):
                if file_path.is_file():
                    arcname = str(file_path.relative_to(tmp_path))
                    zf.write(file_path, arcname)


def repair_wheel_linux(wheel: Path, dest_dir: Path) -> Path:
    """Repair wheel using auditwheel (Linux)."""
    # auditwheel 6.5.0+ automatically generates SBOM for bundled libs
    subprocess.run(
        ["auditwheel", "repair", "-w", str(dest_dir), str(wheel)],
        check=True,
    )
    # Find the repaired wheel
    repaired = list(dest_dir.glob("*.whl"))
    if not repaired:
        raise RuntimeError("No repaired wheel found")
    return repaired[0]


def repair_wheel_windows(wheel: Path, dest_dir: Path) -> Path:
    """Repair wheel using delvewheel (Windows)."""
    subprocess.run(
        ["delvewheel", "repair", "-w", str(dest_dir), str(wheel)],
        check=True,
    )
    # Find the repaired wheel
    repaired = list(dest_dir.glob("*.whl"))
    if not repaired:
        raise RuntimeError("No repaired wheel found")
    return repaired[0]


def repair_wheel_macos(wheel: Path, dest_dir: Path) -> Path:
    """Copy wheel for macOS (no repair needed)."""
    dest = dest_dir / wheel.name
    shutil.copy(wheel, dest)
    return dest


def main() -> None:
    """Main entry point."""
    parser = argparse.ArgumentParser(description="Repair wheel and inject SBOM")
    parser.add_argument("wheel", type=Path, help="Wheel file to repair")
    parser.add_argument("dest_dir", type=Path, help="Destination directory")
    args = parser.parse_args()

    args.dest_dir.mkdir(parents=True, exist_ok=True)

    print(f"Repairing wheel: {args.wheel}")
    print(f"Platform: {sys.platform}")
    print(f"Destination: {args.dest_dir}")

    # Check if this is a Pyodide wheel (built in Linux container but for wasm32)
    wheel_name = args.wheel.name.lower()
    is_pyodide = "pyodide" in wheel_name or "wasm32" in wheel_name

    # Step 1: Run platform-specific repair
    if is_pyodide:
        # Pyodide wheels are already repaired by auditwheel-emscripten
        # Just copy and inject SBOM
        print("Detected Pyodide wheel, skipping native repair")
        repaired_wheel = repair_wheel_macos(args.wheel, args.dest_dir)
    elif sys.platform == "linux":
        repaired_wheel = repair_wheel_linux(args.wheel, args.dest_dir)
    elif sys.platform in ["win32", "cygwin"]:
        repaired_wheel = repair_wheel_windows(args.wheel, args.dest_dir)
    elif sys.platform == "darwin":
        repaired_wheel = repair_wheel_macos(args.wheel, args.dest_dir)
    else:
        raise RuntimeError(f"Unsupported platform: {sys.platform}")

    print(f"Repaired wheel: {repaired_wheel}")

    # Step 2: Generate and inject SBOM
    script_dir = Path(__file__).parent
    sbom_script = script_dir / "generate_sbom.py"

    # Get version from wheel name (e.g., pandas-3.0.0-cp311-...)
    version = repaired_wheel.stem.split("-")[1]

    # Generate SBOM to temp file
    with tempfile.NamedTemporaryFile(suffix=".cdx.json", delete=False) as f:
        sbom_path = Path(f.name)

    final_sbom = sbom_path.parent / "pandas.cdx.json"
    try:
        print(
            f"Running: {sys.executable} {sbom_script} {sbom_path} --version {version}"
        )
        subprocess.run(
            [sys.executable, str(sbom_script), str(sbom_path), "--version", version],
            check=True,
        )
        print(f"Generated SBOM: {sbom_path}")

        sbom_path.rename(final_sbom)

        print(f"Injecting SBOM into {repaired_wheel}")
        inject_sbom(repaired_wheel, final_sbom)
        print("SBOM injection complete")
    finally:
        if sbom_path.exists():
            sbom_path.unlink()
        if final_sbom.exists():
            final_sbom.unlink()


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
Generate a CycloneDX SBOM for pandas vendored components.

This script generates a Software Bill of Materials (SBOM) in CycloneDX 1.6 format
documenting code that pandas has derived from or incorporates from other projects.
This is in compliance with PEP 770.

The vendored components are defined in LICENSES/vendored.toml.

Usage:
    python scripts/generate_sbom.py output_path
    python scripts/generate_sbom.py -  # print to stdout

To validate the generated SBOM:
    check-jsonschema --schemafile \
        https://cyclonedx.org/schema/bom-1.6.schema.json output.json
"""

import argparse
from datetime import (
    datetime,
    timezone,
)
import hashlib
import json
import os
from pathlib import Path
import tomllib


def is_spdx_expression(license_str: str) -> bool:
    """Check if a license string is an SPDX expression (vs a single ID)."""
    # SPDX expressions contain operators like OR, AND, WITH
    return any(op in license_str for op in (" OR ", " AND ", " WITH "))


def load_vendored_components(manifest_path: Path | None = None) -> list[dict]:
    """Load vendored components from LICENSES/vendored.toml manifest."""
    if manifest_path is None:
        # Default to LICENSES/vendored.toml relative to repo root
        repo_root = Path(__file__).parent.parent
        manifest_path = repo_root / "LICENSES" / "vendored.toml"

    with manifest_path.open("rb") as f:
        manifest = tomllib.load(f)

    components = []
    for comp in manifest.get("component", []):
        license_str = comp["license"]
        components.append(
            {
                "name": comp["name"],
                "bom_ref": f"{comp['name']}-derived",
                "description": comp["description"],
                "license": license_str,
                "is_expression": is_spdx_expression(license_str),
                "purl": comp["purl"],
                "website": comp["website"],
            }
        )
    return components


def get_pandas_version() -> str:
    """Get the pandas version from installed package.

    During CI wheel builds, use --version flag instead since pandas
    is not installed at that point.
    """
    try:
        from pandas import __version__

        return __version__
    except ImportError:
        # Return placeholder if pandas is not installed
        return "0.0.0.dev0"


def _reproducible_timestamp() -> str:
    """Return an ISO-8601 timestamp honoring SOURCE_DATE_EPOCH if set.

    Matches the reproducible-builds convention already honored by
    meson-python for wheel file mtimes. Falls back to wall-clock UTC
    only when no SOURCE_DATE_EPOCH is provided.
    """
    sde = os.environ.get("SOURCE_DATE_EPOCH")
    if sde:
        return datetime.fromtimestamp(int(sde), timezone.utc).isoformat()
    return datetime.now(timezone.utc).isoformat()


def _deterministic_serial(version: str, manifest_path: Path) -> str:
    """Build a urn:uuid serialNumber deterministic in pandas version + manifest.

    CycloneDX requires serialNumber to be unique per BOM, but for
    reproducible builds we derive it from inputs so repeated invocations
    produce byte-identical output. Hashing manifest bytes + pandas
    version yields a stable UUID that still changes when either input
    changes.
    """
    manifest_bytes = manifest_path.read_bytes()
    digest = hashlib.sha256(manifest_bytes + version.encode("utf-8")).hexdigest()
    # Lay out the 32-hex digest as a canonical UUID string (8-4-4-4-12).
    u = f"{digest[0:8]}-{digest[8:12]}-{digest[12:16]}-{digest[16:20]}-{digest[20:32]}"
    return f"urn:uuid:{u}"


def generate_sbom(
    version: str | None = None, manifest_path: Path | None = None
) -> dict:
    """Generate the CycloneDX SBOM document."""
    if version is None:
        version = get_pandas_version()
    if manifest_path is None:
        manifest_path = Path(__file__).parent.parent / "LICENSES" / "vendored.toml"

    vendored_components = load_vendored_components(manifest_path)

    timestamp = _reproducible_timestamp()
    serial_number = _deterministic_serial(version, manifest_path)

    # Build components list
    components = []
    dependency_refs = []

    for comp in vendored_components:
        # CycloneDX uses "expression" for SPDX expressions, "id" for single license
        if comp["is_expression"]:
            license_entry = {"expression": comp["license"]}
        else:
            license_entry = {"license": {"id": comp["license"]}}

        component = {
            "type": "library",
            "bom-ref": comp["bom_ref"],
            "name": comp["name"],
            "description": comp["description"],
            "licenses": [license_entry],
            "purl": comp["purl"],
            "externalReferences": [
                {"type": "website", "url": comp["website"]},
            ],
        }
        components.append(component)
        dependency_refs.append(comp["bom_ref"])

    # Single bom-ref shared by metadata.component and dependencies[0]
    # so CycloneDX consumers can resolve the root of the dependency
    # graph. See
    # https://cyclonedx.org/use-cases/software-dependencies/.
    root_bom_ref = f"pkg:pypi/pandas@{version}"

    sbom = {
        "$schema": "https://cyclonedx.org/schema/bom-1.6.schema.json",
        "bomFormat": "CycloneDX",
        "specVersion": "1.6",
        "serialNumber": serial_number,
        "version": 1,
        "metadata": {
            "timestamp": timestamp,
            "tools": {
                "components": [
                    {
                        "type": "application",
                        "name": "pandas-sbom-generator",
                        "version": "1.0.0",
                    }
                ]
            },
            "component": {
                "type": "library",
                "bom-ref": root_bom_ref,
                "name": "pandas",
                "version": version,
                "purl": root_bom_ref,
                "description": "Powerful data structures for data analysis, "
                "time series, and statistics",
                "licenses": [{"license": {"id": "BSD-3-Clause"}}],
                "externalReferences": [
                    {"type": "website", "url": "https://pandas.pydata.org"},
                    {
                        "type": "vcs",
                        "url": "https://github.com/pandas-dev/pandas",
                    },
                    {
                        "type": "documentation",
                        "url": "https://pandas.pydata.org/docs/",
                    },
                ],
            },
        },
        "components": components,
        "dependencies": [
            {
                "ref": root_bom_ref,
                "dependsOn": dependency_refs,
            }
        ],
    }

    return sbom


def main() -> None:
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Generate CycloneDX SBOM for pandas vendored components"
    )
    parser.add_argument(
        "output",
        nargs="?",
        default="-",
        help="Output file path (use '-' for stdout, default: stdout)",
    )
    parser.add_argument(
        "--version",
        help="Override pandas version (default: auto-detect)",
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        help="Path to vendored.toml manifest (default: LICENSES/vendored.toml)",
    )
    args = parser.parse_args()

    sbom = generate_sbom(version=args.version, manifest_path=args.manifest)
    sbom_json = json.dumps(sbom, indent=2, ensure_ascii=False)

    if args.output == "-":
        print(sbom_json)
    else:
        output_path = Path(args.output)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(sbom_json, encoding="utf-8")


if __name__ == "__main__":
    main()

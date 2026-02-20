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
import json
from pathlib import Path
import tomllib
import uuid


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


def generate_sbom(
    version: str | None = None, manifest_path: Path | None = None
) -> dict:
    """Generate the CycloneDX SBOM document."""
    if version is None:
        version = get_pandas_version()

    vendored_components = load_vendored_components(manifest_path)

    timestamp = datetime.now(timezone.utc).isoformat()
    serial_number = f"urn:uuid:{uuid.uuid4()}"

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
                "name": "pandas",
                "version": version,
                "purl": f"pkg:pypi/pandas@{version}",
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
                "ref": f"pkg:pypi/pandas@{version}",
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

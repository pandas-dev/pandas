# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.

"""OpenAPI spec utils."""
from __future__ import annotations

import os
import typing
from pathlib import Path

if typing.TYPE_CHECKING:
    from openapi_core.spec.paths import Spec

HERE = Path(os.path.dirname(__file__)).resolve()


def get_openapi_spec() -> Spec:
    """Get the OpenAPI spec object."""
    from openapi_core.spec.paths import Spec

    openapi_spec_dict = get_openapi_spec_dict()
    return Spec.from_dict(openapi_spec_dict)  # type:ignore[arg-type]


def get_openapi_spec_dict() -> dict[str, typing.Any]:
    """Get the OpenAPI spec as a dictionary."""
    from ruamel.yaml import YAML

    path = HERE / "rest-api.yml"
    yaml = YAML(typ="safe")
    return yaml.load(path.read_text(encoding="utf-8"))

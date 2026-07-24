"""Notebook format validators."""

# Copyright (c) IPython Development Team.
# Distributed under the terms of the Modified BSD License.
from __future__ import annotations

import json
import pprint
import warnings
from copy import deepcopy
from pathlib import Path
from textwrap import dedent
from typing import Any, Optional

from ._imports import import_item
from .corpus.words import generate_corpus_id
from .json_compat import ValidationError, _validator_for_name, get_current_validator
from .reader import get_version
from .warnings import DuplicateCellId, MissingIDFieldWarning

validators = {}
_deprecated = object()


__all__ = [
    "ValidationError",
    "get_validator",
    "isvalid",
    "NotebookValidationError",
    "better_validation_error",
    "normalize",
    "validate",
    "iter_validate",
]


def _relax_additional_properties(obj):
    """relax any `additionalProperties`"""
    if isinstance(obj, dict):
        for key, value in obj.items():
            value = (  # noqa: PLW2901
                True if key == "additionalProperties" else _relax_additional_properties(value)
            )
            obj[key] = value
    elif isinstance(obj, list):
        for i, value in enumerate(obj):
            obj[i] = _relax_additional_properties(value)
    return obj


def _allow_undefined(schema):
    schema["definitions"]["cell"]["oneOf"].append({"$ref": "#/definitions/unrecognized_cell"})
    schema["definitions"]["output"]["oneOf"].append({"$ref": "#/definitions/unrecognized_output"})
    return schema


def get_validator(version=None, version_minor=None, relax_add_props=False, name=None):
    """Load the JSON schema into a Validator"""
    if version is None:
        from . import current_nbformat

        version = current_nbformat

    v = import_item("nbformat.v%s" % version)
    current_minor = getattr(v, "nbformat_minor", 0)
    if version_minor is None:
        version_minor = current_minor

    current_validator = _validator_for_name(name) if name else get_current_validator()

    version_tuple = (current_validator.name, version, version_minor)

    if version_tuple not in validators:
        try:
            schema_json = _get_schema_json(v, version=version, version_minor=version_minor)
        except AttributeError:
            return None

        if current_minor < version_minor:
            # notebook from the future, relax all `additionalProperties: False` requirements
            schema_json = _relax_additional_properties(schema_json)
            # and allow undefined cell types and outputs
            schema_json = _allow_undefined(schema_json)

        validators[version_tuple] = current_validator(schema_json)

    if relax_add_props:
        try:
            schema_json = _get_schema_json(v, version=version, version_minor=version_minor)
        except AttributeError:
            return None

        # this allows properties to be added for intermediate
        # representations while validating for all other kinds of errors
        schema_json = _relax_additional_properties(schema_json)
        validators[version_tuple] = current_validator(schema_json)

    return validators[version_tuple]


def _get_schema_json(v, version=None, version_minor=None):
    """
    Gets the json schema from a given imported library and nbformat version.
    """
    if (version, version_minor) in v.nbformat_schema:
        schema_path = str(Path(v.__file__).parent / v.nbformat_schema[(version, version_minor)])
    elif version_minor > v.nbformat_minor:
        # load the latest schema
        schema_path = str(Path(v.__file__).parent / v.nbformat_schema[(None, None)])
    else:
        msg = "Cannot find appropriate nbformat schema file."
        raise AttributeError(msg)
    with Path(schema_path).open(encoding="utf8") as f:
        schema_json = json.load(f)
    return schema_json  # noqa: RET504


def isvalid(nbjson, ref=None, version=None, version_minor=None):
    """Checks whether the given notebook JSON conforms to the current
    notebook format schema. Returns True if the JSON is valid, and
    False otherwise.

    To see the individual errors that were encountered, please use the
    `validate` function instead.
    """
    orig = deepcopy(nbjson)
    try:
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=DeprecationWarning)
            warnings.filterwarnings("ignore", category=MissingIDFieldWarning)
            validate(nbjson, ref, version, version_minor, repair_duplicate_cell_ids=False)
    except ValidationError:
        return False
    else:
        return True
    finally:
        if nbjson != orig:
            raise AssertionError


def _format_as_index(indices):
    """
    (from jsonschema._utils.format_as_index, copied to avoid relying on private API)

    Construct a single string containing indexing operations for the indices.

    For example, [1, 2, "foo"] -> [1][2]["foo"]
    """

    if not indices:
        return ""
    return "[%s]" % "][".join(repr(index) for index in indices)


_ITEM_LIMIT = 16
_STR_LIMIT = 64


def _truncate_obj(obj):
    """Truncate objects for use in validation tracebacks

    Cell and output lists are squashed, as are long strings, lists, and dicts.
    """
    if isinstance(obj, dict):
        truncated_dict = {k: _truncate_obj(v) for k, v in list(obj.items())[:_ITEM_LIMIT]}
        if isinstance(truncated_dict.get("cells"), list):
            truncated_dict["cells"] = ["...%i cells..." % len(obj["cells"])]
        if isinstance(truncated_dict.get("outputs"), list):
            truncated_dict["outputs"] = ["...%i outputs..." % len(obj["outputs"])]

        if len(obj) > _ITEM_LIMIT:
            truncated_dict["..."] = "%i keys truncated" % (len(obj) - _ITEM_LIMIT)
        return truncated_dict
    if isinstance(obj, list):
        truncated_list = [_truncate_obj(item) for item in obj[:_ITEM_LIMIT]]
        if len(obj) > _ITEM_LIMIT:
            truncated_list.append("...%i items truncated..." % (len(obj) - _ITEM_LIMIT))
        return truncated_list
    if isinstance(obj, str):
        truncated_str = obj[:_STR_LIMIT]
        if len(obj) > _STR_LIMIT:
            truncated_str += "..."
        return truncated_str
    return obj


class NotebookValidationError(ValidationError):  # type:ignore[misc]
    """Schema ValidationError with truncated representation

    to avoid massive verbose tracebacks.
    """

    def __init__(self, original, ref=None):
        """Initialize the error class."""
        self.original = original
        self.ref = getattr(self.original, "ref", ref)
        self.message = self.original.message

    def __getattr__(self, key):
        """Get an attribute from the error."""
        return getattr(self.original, key)

    def __unicode__(self):
        """Custom str for validation errors

        avoids dumping full schema and notebook to logs
        """
        error = self.original
        instance = _truncate_obj(error.instance)

        return "\n".join(
            [
                error.message,
                "",
                "Failed validating {!r} in {}{}:".format(
                    error.validator,
                    self.ref or "notebook",
                    _format_as_index(list(error.relative_schema_path)[:-1]),
                ),
                "",
                "On instance%s:" % _format_as_index(error.relative_path),
                pprint.pformat(instance, width=78),
            ]
        )

    __str__ = __unicode__


def better_validation_error(error, version, version_minor):
    """Get better ValidationError on oneOf failures

    oneOf errors aren't informative.
    if it's a cell type or output_type error,
    try validating directly based on the type for a better error message
    """
    if not len(error.schema_path):
        return error
    key = error.schema_path[-1]
    ref = None
    if key.endswith("Of"):
        if isinstance(error.instance, dict):
            if "cell_type" in error.instance:
                ref = error.instance["cell_type"] + "_cell"
            elif "output_type" in error.instance:
                ref = error.instance["output_type"]

        if ref:
            try:
                validate(
                    error.instance,
                    ref,
                    version=version,
                    version_minor=version_minor,
                )
            except ValidationError as sub_error:
                # keep extending relative path
                error.relative_path.extend(sub_error.relative_path)
                sub_error.relative_path = error.relative_path
                better = better_validation_error(sub_error, version, version_minor)
                if better.ref is None:
                    better.ref = ref
                return better
            except Exception:  # noqa: S110
                # if it fails for some reason,
                # let the original error through
                pass
    return NotebookValidationError(error, ref)


def normalize(
    nbdict: Any,
    version: Optional[int] = None,
    version_minor: Optional[int] = None,
    *,
    relax_add_props: bool = False,
    strip_invalid_metadata: bool = False,
) -> tuple[int, Any]:
    """
    Normalise a notebook prior to validation.

    This tries to implement a couple of normalisation steps to standardise
    notebooks and make validation easier.

    You should in general not rely on this function and make sure the notebooks
    that reach nbformat are already in a normal form. If not you likely have a bug,
    and may have security issues.

    Parameters
    ----------
    nbdict : dict
        notebook document
    version : int
    version_minor : int
    relax_add_props : bool
        Whether to allow extra property in the Json schema validating the
        notebook.
    strip_invalid_metadata : bool
        Whether to strip metadata that does not exist in the Json schema when
        validating the notebook.

    Returns
    -------
    changes : int
        number of changes in the notebooks
    notebook : dict
        deep-copy of the original object with relevant changes.

    """
    nbdict = deepcopy(nbdict)
    nbdict_version, nbdict_version_minor = get_version(nbdict)
    if version is None:
        version = nbdict_version
    if version_minor is None:
        version_minor = nbdict_version_minor
    return _normalize(
        nbdict,
        version,
        version_minor,
        True,
        relax_add_props=relax_add_props,
        strip_invalid_metadata=strip_invalid_metadata,
    )


def _normalize(
    nbdict: Any,
    version: int,
    version_minor: int,
    repair_duplicate_cell_ids: bool,
    relax_add_props: bool,
    strip_invalid_metadata: bool,
) -> tuple[int, Any]:
    """
    Private normalisation routine.

    This function attempts to normalize the `nbdict` passed to it.

    As `_normalize()` is currently used both in `validate()` (for
    historical reasons), and in the `normalize()` public function,
    `_normalize()` does currently mutate `nbdict`.
    Ideally, once `validate()` stops calling `_normalize()`, `_normalize()`
    may stop mutating `nbdict`.

    """
    changes = 0

    if (version, version_minor) >= (4, 5):
        # if we support cell ids ensure default ids are provided
        for cell in nbdict["cells"]:
            if "id" not in cell:
                warnings.warn(
                    "Cell is missing an id field, this will become"
                    " a hard error in future nbformat versions. You may want"
                    " to use `normalize()` on your notebooks before validations"
                    " (available since nbformat 5.1.4). Previous versions of nbformat"
                    " are fixing this issue transparently, and will stop doing so"
                    " in the future.",
                    MissingIDFieldWarning,
                    stacklevel=3,
                )
                # Generate cell ids if any are missing
                if repair_duplicate_cell_ids:
                    cell["id"] = generate_corpus_id()
                    changes += 1

        # if we support cell ids check for uniqueness when validating the whole notebook
        seen_ids = set()
        for cell in nbdict["cells"]:
            if "id" not in cell:
                continue
            cell_id = cell["id"]
            if cell_id in seen_ids:
                # Best effort to repair if we find a duplicate id
                if repair_duplicate_cell_ids:
                    new_id = generate_corpus_id()
                    cell["id"] = new_id
                    changes += 1
                    warnings.warn(
                        f"Non-unique cell id {cell_id!r} detected. Corrected to {new_id!r}.",
                        DuplicateCellId,
                        stacklevel=3,
                    )
                else:
                    msg = f"Non-unique cell id '{cell_id}' detected."
                    raise ValidationError(msg)
            seen_ids.add(cell_id)
    if strip_invalid_metadata:
        changes += _strip_invalida_metadata(
            nbdict, version, version_minor, relax_add_props=relax_add_props
        )
    return changes, nbdict


def _dep_warn(field):
    warnings.warn(
        dedent(
            f"""`{field}` kwargs of validate has been deprecated for security
        reasons, and will be removed soon.

        Please explicitly use the `n_changes, new_notebook = nbformat.validator.normalize(old_notebook, ...)` if you wish to
        normalise your notebook. `normalize` is available since nbformat 5.5.0

        """
        ),
        DeprecationWarning,
        stacklevel=3,
    )


def validate(
    nbdict: Any = None,
    ref: Optional[str] = None,
    version: Optional[int] = None,
    version_minor: Optional[int] = None,
    relax_add_props: bool = False,
    nbjson: Any = None,
    repair_duplicate_cell_ids: bool = _deprecated,  # type: ignore[assignment]
    strip_invalid_metadata: bool = _deprecated,  # type: ignore[assignment]
) -> None:
    """Checks whether the given notebook dict-like object
    conforms to the relevant notebook format schema.

    Parameters
    ----------
    nbdict : dict
        notebook document
    ref : optional, str
        reference to the subset of the schema we want to validate against.
        for example ``"markdown_cell"``, `"code_cell"` ....
    version : int
    version_minor : int
    relax_add_props : bool
        Whether to allow extra properties in the JSON schema validating the notebook.
        When True, all known fields are validated, but unknown fields are ignored.
    nbjson
    repair_duplicate_cell_ids : bool
        Deprecated since 5.5.0 - will be removed in the future.
    strip_invalid_metadata : bool
        Deprecated since 5.5.0 - will be removed in the future.

    Returns
    -------
    None

    Raises
    ------
    ValidationError if not valid.

    Notes
    -----
    Prior to Nbformat 5.5.0 the `validate` and `isvalid` method would silently
    try to fix invalid notebook and mutate arguments. This behavior is deprecated
    and will be removed in a near future.

    Please explicitly call `normalize` if you need to normalize notebooks.
    """
    assert isinstance(ref, str) or ref is None

    if strip_invalid_metadata is _deprecated:
        strip_invalid_metadata = False
    else:
        _dep_warn("strip_invalid_metadata")

    if repair_duplicate_cell_ids is _deprecated:
        repair_duplicate_cell_ids = True
    else:
        _dep_warn("repair_duplicate_cell_ids")

    # backwards compatibility for nbjson argument
    if nbdict is not None:
        pass
    elif nbjson is not None:
        nbdict = nbjson
    else:
        msg = "validate() missing 1 required argument: 'nbdict'"
        raise TypeError(msg)

    if ref is None:
        # if ref is not specified, we have a whole notebook, so we can get the version
        nbdict_version, nbdict_version_minor = get_version(nbdict)
        if version is None:
            version = nbdict_version
        if version_minor is None:
            version_minor = nbdict_version_minor
    # if ref is specified, and we don't have a version number, assume we're validating against 1.0
    elif version is None:
        version, version_minor = 1, 0

    if ref is None:
        assert isinstance(version, int)
        assert isinstance(version_minor, int)
        _normalize(
            nbdict,
            version,
            version_minor,
            repair_duplicate_cell_ids,
            relax_add_props=relax_add_props,
            strip_invalid_metadata=strip_invalid_metadata,
        )

    for error in iter_validate(
        nbdict,
        ref=ref,
        version=version,
        version_minor=version_minor,
        relax_add_props=relax_add_props,
        strip_invalid_metadata=strip_invalid_metadata,
    ):
        raise error


def _get_errors(
    nbdict: Any, version: int, version_minor: int, relax_add_props: bool, *args: Any
) -> Any:
    validator = get_validator(version, version_minor, relax_add_props=relax_add_props)
    if not validator:
        msg = f"No schema for validating v{version}.{version_minor} notebooks"
        raise ValidationError(msg)
    iter_errors = validator.iter_errors(nbdict, *args)
    errors = list(iter_errors)
    # jsonschema gives the best error messages.
    if len(errors) and validator.name != "jsonschema":
        validator = get_validator(
            version=version,
            version_minor=version_minor,
            relax_add_props=relax_add_props,
            name="jsonschema",
        )
        return validator.iter_errors(nbdict, *args)
    return iter(errors)


def _strip_invalida_metadata(
    nbdict: Any, version: int, version_minor: int, relax_add_props: bool
) -> int:
    """
    This function tries to extract metadata errors from the validator and fix
    them if necessary. This mostly mean stripping unknown keys from metadata
    fields, or removing metadata fields altogether.

    Parameters
    ----------
    nbdict : dict
        notebook document
    version : int
    version_minor : int
    relax_add_props : bool
        Whether to allow extra property in the Json schema validating the
        notebook.

    Returns
    -------
    int
        number of modifications

    """
    errors = _get_errors(nbdict, version, version_minor, relax_add_props)
    changes = 0
    if len(list(errors)) > 0:
        # jsonschema gives a better error tree.
        validator = get_validator(
            version=version,
            version_minor=version_minor,
            relax_add_props=relax_add_props,
            name="jsonschema",
        )
        if not validator:
            msg = f"No jsonschema for validating v{version}.{version_minor} notebooks"
            raise ValidationError(msg)
        errors = validator.iter_errors(nbdict)
        error_tree = validator.error_tree(errors)
        if "metadata" in error_tree:
            for key in error_tree["metadata"]:
                nbdict["metadata"].pop(key, None)
                changes += 1

        if "cells" in error_tree:
            number_of_cells = len(nbdict.get("cells", 0))
            for cell_idx in range(number_of_cells):
                # Cells don't report individual metadata keys as having failed validation
                # Instead it reports that it failed to validate against each cell-type definition.
                # We have to delve into why those definitions failed to uncover which metadata
                # keys are misbehaving.
                if "oneOf" in error_tree["cells"][cell_idx].errors:
                    intended_cell_type = nbdict["cells"][cell_idx]["cell_type"]
                    schemas_by_index = [
                        ref["$ref"]
                        for ref in error_tree["cells"][cell_idx].errors["oneOf"].schema["oneOf"]
                    ]
                    cell_type_definition_name = f"#/definitions/{intended_cell_type}_cell"
                    if cell_type_definition_name in schemas_by_index:
                        schema_index = schemas_by_index.index(cell_type_definition_name)
                        for error in error_tree["cells"][cell_idx].errors["oneOf"].context:
                            rel_path = error.relative_path
                            error_for_intended_schema = error.schema_path[0] == schema_index
                            is_top_level_metadata_key = (
                                len(rel_path) == 2 and rel_path[0] == "metadata"
                            )
                            if error_for_intended_schema and is_top_level_metadata_key:
                                nbdict["cells"][cell_idx]["metadata"].pop(rel_path[1], None)
                                changes += 1

    return changes


def iter_validate(
    nbdict=None,
    ref=None,
    version=None,
    version_minor=None,
    relax_add_props=False,
    nbjson=None,
    strip_invalid_metadata=False,
):
    """Checks whether the given notebook dict-like object conforms to the
    relevant notebook format schema.

    Returns a generator of all ValidationErrors if not valid.

    Notes
    -----
    To fix: For security reasons, this function should *never* mutate its `nbdict` argument, and
    should *never* try to validate a mutated or modified version of its notebook.

    """
    # backwards compatibility for nbjson argument
    if nbdict is not None:
        pass
    elif nbjson is not None:
        nbdict = nbjson
    else:
        msg = "iter_validate() missing 1 required argument: 'nbdict'"
        raise TypeError(msg)

    if version is None:
        version, version_minor = get_version(nbdict)

    if ref:
        try:
            errors = _get_errors(
                nbdict,
                version,
                version_minor,
                relax_add_props,
                {"$ref": "#/definitions/%s" % ref},
            )
        except ValidationError as e:
            yield e
            return

    else:
        if strip_invalid_metadata:
            _strip_invalida_metadata(nbdict, version, version_minor, relax_add_props)

        # Validate one more time to ensure that us removing metadata
        # didn't cause another complex validation issue in the schema.
        # Also to ensure that higher-level errors produced by individual metadata validation
        # failures are removed.
        try:
            errors = _get_errors(nbdict, version, version_minor, relax_add_props)
        except ValidationError as e:
            yield e
            return

    for error in errors:
        yield better_validation_error(error, version, version_minor)

"""Extension to create automatic documentation from code docstrings.

Automatically insert docstrings for functions, classes or whole modules into
the doctree, thus avoiding duplication between docstrings and documentation
for those who like elaborate docstrings.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import sphinx
from sphinx.config import ENUM
from sphinx.ext.autodoc._directive import AutodocDirective
from sphinx.ext.autodoc._dynamic._importer import (
    _import_module as import_module,
)
from sphinx.ext.autodoc._dynamic._mock import ismock, mock, undecorate
from sphinx.ext.autodoc._event_listeners import between, cut_lines
from sphinx.ext.autodoc._legacy_class_based._directive_options import (
    Options,
    annotation_option,
    bool_option,
    class_doc_from_option,
    exclude_members_option,
    identity,
    inherited_members_option,
    member_order_option,
    members_option,
    merge_members_option,
)
from sphinx.ext.autodoc._legacy_class_based._documenters import (
    AttributeDocumenter,
    ClassDocumenter,
    ClassLevelDocumenter,
    DataDocumenter,
    DataDocumenterMixinBase,
    DecoratorDocumenter,
    DocstringSignatureMixin,
    DocstringStripSignatureMixin,
    Documenter,
    ExceptionDocumenter,
    FunctionDocumenter,
    GenericAliasMixin,
    MethodDocumenter,
    ModuleDocumenter,
    ModuleLevelDocumenter,
    NonDataDescriptorMixin,
    ObjectMember,
    PropertyDocumenter,
    RuntimeInstanceAttributeMixin,
    SlotsMixin,
    UninitializedGlobalVariableMixin,
    UninitializedInstanceAttributeMixin,
    autodoc_attrgetter,
    py_ext_sig_re,
    special_member_re,
)
from sphinx.ext.autodoc._legacy_class_based._sentinels import (
    ALL,
    EMPTY,
    INSTANCEATTR,
    SLOTSATTR,
    SUPPRESS,
    UNINITIALIZED_ATTR,
)
from sphinx.ext.autodoc.importer import get_class_members, import_object
from sphinx.ext.autodoc.typehints import _merge_typehints

if TYPE_CHECKING:
    from sphinx.application import Sphinx
    from sphinx.config import Config
    from sphinx.ext.autodoc._property_types import _AutodocObjType
    from sphinx.util.typing import ExtensionMetadata

__all__ = (
    # Useful event listener factories for autodoc-process-docstring
    'cut_lines',
    'between',
    # Documenters
    'AttributeDocumenter',
    'ClassDocumenter',
    'DataDocumenter',
    'DecoratorDocumenter',
    'ExceptionDocumenter',
    'FunctionDocumenter',
    'MethodDocumenter',
    'ModuleDocumenter',
    'PropertyDocumenter',
    # This class is only used in ``sphinx.ext.autodoc.directive``,
    # but we export it here for compatibility.
    # See: https://github.com/sphinx-doc/sphinx/issues/4538
    'Options',
    # Option spec functions.
    # Exported for compatibility.
    'annotation_option',
    'bool_option',
    'class_doc_from_option',
    'exclude_members_option',
    'identity',
    'inherited_members_option',
    'member_order_option',
    'members_option',
    'merge_members_option',
    # Sentinels.
    # Exported for compatibility.
    'ALL',
    'EMPTY',
    'INSTANCEATTR',
    'SLOTSATTR',
    'SUPPRESS',
    'UNINITIALIZED_ATTR',
    # Miscellaneous other names.
    # Exported for compatibility.
    'ObjectMember',
    'py_ext_sig_re',
    'special_member_re',
    'ModuleLevelDocumenter',
    'ClassLevelDocumenter',
    'DocstringSignatureMixin',
    'DocstringStripSignatureMixin',
    'DataDocumenterMixinBase',
    'GenericAliasMixin',
    'UninitializedGlobalVariableMixin',
    'NonDataDescriptorMixin',
    'SlotsMixin',
    'RuntimeInstanceAttributeMixin',
    'UninitializedInstanceAttributeMixin',
    'autodoc_attrgetter',
    'Documenter',
    'get_class_members',
    'import_module',
    'import_object',
    'ismock',
    'mock',
    'undecorate',
)


def setup(app: Sphinx) -> ExtensionMetadata:
    app.add_config_value(
        'autoclass_content',
        'class',
        'env',
        types=ENUM('both', 'class', 'init'),
    )
    app.add_config_value(
        'autodoc_member_order',
        'alphabetical',
        'env',
        types=ENUM('alphabetical', 'bysource', 'groupwise'),
    )
    app.add_config_value(
        'autodoc_class_signature',
        'mixed',
        'env',
        types=ENUM('mixed', 'separated'),
    )
    app.add_config_value('autodoc_default_options', {}, 'env', types=frozenset({dict}))
    app.add_config_value(
        'autodoc_docstring_signature', True, 'env', types=frozenset({bool})
    )
    app.add_config_value(
        'autodoc_mock_imports', [], 'env', types=frozenset({list, tuple})
    )
    app.add_config_value(
        'autodoc_typehints',
        'signature',
        'env',
        types=ENUM('signature', 'description', 'none', 'both'),
    )
    app.add_config_value(
        'autodoc_typehints_description_target',
        'all',
        'env',
        types=ENUM('all', 'documented', 'documented_params'),
    )
    app.add_config_value('autodoc_type_aliases', {}, 'env', types=frozenset({dict}))
    app.add_config_value(
        'autodoc_typehints_format',
        'short',
        'env',
        types=ENUM('fully-qualified', 'short'),
    )
    app.add_config_value('autodoc_warningiserror', True, 'env', types=frozenset({bool}))
    app.add_config_value(
        'autodoc_inherit_docstrings', True, 'env', types=frozenset({bool})
    )
    app.add_config_value(
        'autodoc_preserve_defaults', False, 'env', types=frozenset({bool})
    )
    app.add_config_value(
        'autodoc_use_type_comments', True, 'env', types=frozenset({bool})
    )
    app.add_config_value(
        'autodoc_use_legacy_class_based', False, 'env', types=frozenset({bool})
    )

    app.add_event('autodoc-before-process-signature')
    app.add_event('autodoc-process-docstring')
    app.add_event('autodoc-process-signature')
    app.add_event('autodoc-skip-member')
    app.add_event('autodoc-process-bases')

    app.connect('object-description-transform', _merge_typehints)

    app.connect('config-inited', _register_directives)

    return {
        'version': sphinx.__display_version__,
        'parallel_read_safe': True,
    }


def _register_directives(app: Sphinx, config: Config) -> None:
    if not config.autodoc_use_legacy_class_based:
        obj_type: _AutodocObjType
        for obj_type in (
            'module',
            'class',
            'exception',
            'function',
            'decorator',
            'method',
            'property',
            'attribute',
            'data',
            'type',
        ):
            # register the automodule, autoclass, etc. directives
            app.add_directive(f'auto{obj_type}', AutodocDirective)
    else:
        from sphinx.ext.autodoc.preserve_defaults import update_defvalue
        from sphinx.ext.autodoc.type_comment import (
            update_annotations_using_type_comments,
        )
        from sphinx.ext.autodoc.typehints import record_typehints

        app.add_autodocumenter(ModuleDocumenter)
        app.add_autodocumenter(ClassDocumenter)
        app.add_autodocumenter(ExceptionDocumenter)
        app.add_autodocumenter(DataDocumenter)
        app.add_autodocumenter(FunctionDocumenter)
        app.add_autodocumenter(DecoratorDocumenter)
        app.add_autodocumenter(MethodDocumenter)
        app.add_autodocumenter(AttributeDocumenter)
        app.add_autodocumenter(PropertyDocumenter)

        app.connect('autodoc-before-process-signature', update_defvalue)
        app.connect(
            'autodoc-before-process-signature', update_annotations_using_type_comments
        )
        app.connect('autodoc-process-signature', record_typehints)

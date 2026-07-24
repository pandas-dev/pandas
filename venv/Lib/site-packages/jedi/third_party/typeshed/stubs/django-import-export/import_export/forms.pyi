from collections.abc import Iterable, Sequence
from typing import Any

from django import forms

from .formats.base_formats import Format
from .resources import ModelResource, Resource

class ImportExportFormBase(forms.Form):
    resource: forms.ChoiceField
    format: forms.ChoiceField
    def __init__(
        self, formats: list[type[Format]], resources: list[type[Resource[Any]]] | None = None, **kwargs: Any
    ) -> None: ...

class ImportForm(ImportExportFormBase):
    import_file: forms.FileField
    field_order: Sequence[str]
    @property
    def media(self) -> forms.Media: ...

class ConfirmImportForm(forms.Form):
    import_file_name: forms.CharField
    original_file_name: forms.CharField
    resource: forms.CharField
    def clean_import_file_name(self) -> str: ...

class ExportForm(ImportExportFormBase):
    export_items: forms.MultipleChoiceField

class SelectableFieldsExportForm(ExportForm):
    resources: Iterable[ModelResource[Any]]
    is_selectable_fields_form: bool
    resource_fields: dict[str, list[str]]
    @staticmethod
    def create_boolean_field_name(resource: ModelResource[Any], field_name: str) -> str: ...
    def get_selected_resource(self) -> ModelResource[Any]: ...
    def get_selected_resource_export_fields(self) -> list[str]: ...

from collections.abc import Sequence
from typing import Any, Generic, TypeVar

from django.db.models import Model

from .instance_loaders import BaseInstanceLoader

_ModelT = TypeVar("_ModelT", bound=Model)

class ResourceOptions(Generic[_ModelT]):
    model: _ModelT | str
    fields: Sequence[str] | None
    exclude: Sequence[str] | None
    instance_loader_class: type[BaseInstanceLoader] | None
    import_id_fields: Sequence[str]
    export_order: Sequence[str] | None
    import_order: Sequence[str] | None
    widgets: dict[str, Any] | None
    use_transactions: bool | None
    skip_unchanged: bool
    report_skipped: bool
    clean_model_instances: bool
    chunk_size: int | None
    skip_diff: bool
    skip_html_diff: bool
    use_bulk: bool
    batch_size: int
    force_init_instance: bool
    using_db: str | None
    store_row_values: bool
    store_instance: bool
    use_natural_foreign_keys: bool

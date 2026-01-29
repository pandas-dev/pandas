from typing import (
    Any,
)

from pandas.compat._optional import import_optional_dependency
from pandas.util._decorators import set_module

from pandas import DataFrame


@set_module("pandas")
def read_iceberg(
    table_identifier: str,
    catalog_name: str | None = None,
    *,
    catalog_properties: dict[str, Any] | None = None,
    columns: list[str] | None = None,
    row_filter: str | None = None,
    case_sensitive: bool = True,
    snapshot_id: int | None = None,
    limit: int | None = None,
    scan_properties: dict[str, Any] | None = None,
) -> DataFrame:
    """
    Read an Apache Iceberg table into a pandas DataFrame.

    .. versionadded:: 3.0.0

    .. warning::

       read_iceberg is experimental and may change without warning.

    Parameters
    ----------
    table_identifier : str
        Table identifier.
    catalog_name : str, optional
        The name of the catalog.
    catalog_properties : dict of {str: str}, optional
        The properties that are used next to the catalog configuration.
    columns : list of str, optional
        A list of strings representing the column names to return in the output
        dataframe.
    row_filter : str, optional
        A string that describes the desired rows.
    case_sensitive : bool, default True
        If True column matching is case sensitive.
    snapshot_id : int, optional
        Snapshot ID to time travel to. By default the table will be scanned as of the
        current snapshot ID.
    limit : int, optional
        An integer representing the number of rows to return in the scan result.
        By default all matching rows will be fetched.
    scan_properties : dict of {str: obj}, optional
        Additional Table properties as a dictionary of string key value pairs to use
        for this scan.

    Returns
    -------
    DataFrame
        DataFrame based on the Iceberg table.

    See Also
    --------
    read_parquet : Read a Parquet file.

    Examples
    --------
    >>> df = pd.read_iceberg(
    ...     table_identifier="my_table",
    ...     catalog_name="my_catalog",
    ...     catalog_properties={"s3.secret-access-key": "my-secret"},
    ...     row_filter="trip_distance >= 10.0",
    ...     columns=["VendorID", "tpep_pickup_datetime"],
    ... )  # doctest: +SKIP
    """
    pyiceberg_catalog = import_optional_dependency("pyiceberg.catalog")
    pyiceberg_expressions = import_optional_dependency("pyiceberg.expressions")
    if catalog_properties is None:
        catalog_properties = {}
    catalog = pyiceberg_catalog.load_catalog(catalog_name, **catalog_properties)
    table = catalog.load_table(table_identifier)
    if row_filter is None:
        row_filter = pyiceberg_expressions.AlwaysTrue()
    if columns is None:
        selected_fields = ("*",)
    else:
        selected_fields = tuple(columns)  # type: ignore[assignment]
    if scan_properties is None:
        scan_properties = {}
    result = table.scan(
        row_filter=row_filter,
        selected_fields=selected_fields,
        case_sensitive=case_sensitive,
        snapshot_id=snapshot_id,
        options=scan_properties,
        limit=limit,
    )
    return result.to_pandas()


def to_iceberg(
    df: DataFrame,
    table_identifier: str,
    catalog_name: str | None = None,
    *,
    catalog_properties: dict[str, Any] | None = None,
    location: str | None = None,
    append: bool = False,
    snapshot_properties: dict[str, str] | None = None,
) -> None:
    """
    Write a DataFrame to an Apache Iceberg table.

    .. versionadded:: 3.0.0

    Parameters
    ----------
    table_identifier : str
        Table identifier.
    catalog_name : str, optional
        The name of the catalog.
    catalog_properties : dict of {str: str}, optional
        The properties that are used next to the catalog configuration.
    location : str, optional
        Location for the table.
    append : bool, default False
        If ``True``, append data to the table, instead of replacing the content.
    snapshot_properties : dict of {str: str}, optional
        Custom properties to be added to the snapshot summary

    See Also
    --------
    read_iceberg : Read an Apache Iceberg table.
    DataFrame.to_parquet : Write a DataFrame in Parquet format.
    """
    pa = import_optional_dependency("pyarrow")
    pyiceberg_catalog = import_optional_dependency("pyiceberg.catalog")
    if catalog_properties is None:
        catalog_properties = {}
    catalog = pyiceberg_catalog.load_catalog(catalog_name, **catalog_properties)
    arrow_table = pa.Table.from_pandas(df)
    table = catalog.create_table_if_not_exists(
        identifier=table_identifier,
        schema=arrow_table.schema,
        location=location,
        # we could add `partition_spec`, `sort_order` and `properties` in the
        # future, but it may not be trivial without exposing PyIceberg objects
    )
    if snapshot_properties is None:
        snapshot_properties = {}
    if append:
        table.append(arrow_table, snapshot_properties=snapshot_properties)
    else:
        table.overwrite(arrow_table, snapshot_properties=snapshot_properties)

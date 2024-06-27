import copy
import datetime
import re
from collections import OrderedDict, namedtuple
from enum import Enum
from typing import Any, Dict, List, Optional, Tuple

from botocore.utils import merge_dicts

SECONDS_IN_ONE_DAY = 24 * 60 * 60
FilterDef = namedtuple(
    "FilterDef",
    [
        # A list of object attributes to check against the filter values.
        # Set to None if filter is not yet implemented in `moto`.
        "attrs_to_check",
        # Description of the filter, e.g. 'Object Identifiers'.
        # Used in filter error messaging.
        "description",
    ],
)


class DbInstanceEngine(str, Enum):
    # https://boto3.amazonaws.com/v1/documentation/api/latest/reference/services/rds/client/create_db_instance.html
    # 2023-11-08
    AURORA_MYSQL = "aurora-mysql"
    AURORA_POSTGRESQL = "aurora-postgresql"
    CUSTOM_ORACLE_EE = "custom-oracle-ee"
    CUSTOM_ORACLE_EE_CDB = "custom-oracle-ee-cdb"
    CUSTOM_SQLSERVER_EE = "custom-sqlserver-ee"
    CUSTOM_SQLSERVER_SE = "custom-sqlserver-se"
    CUSTOM_SQLSERVER_WEB = "custom-sqlserver-web"
    MARIADB = "mariadb"
    MYSQL = "mysql"
    ORACLE_EE = "oracle-ee"
    ORACLE_EE_CDB = "oracle-ee-cdb"
    ORACLE_SE2 = "oracle-se2"
    ORACLE_SE2_CDB = "oracle-se2-cdb"
    POSTGRES = "postgres"
    SQLSERVER_EE = "sqlserver-ee"
    SQLSERVER_SE = "sqlserver-se"
    SQLSERVER_EX = "sqlserver-ex"
    SQLSERVER_WEB = "sqlserver-web"

    @classmethod
    def valid_db_instance_engine(self) -> List[str]:
        return sorted([item.value for item in DbInstanceEngine])


class ClusterEngine(str, Enum):
    AURORA_POSTGRESQL = "aurora-postgresql"
    AURORA_MYSQL = "aurora-mysql"
    RDS_POSTGRESQL = "postgres"
    RDS_MYSQL = "mysql"

    @classmethod
    def list_cluster_engines(self) -> List[str]:
        return sorted([item.value for item in ClusterEngine])

    @classmethod
    def serverless_engines(self) -> List[str]:
        return [ClusterEngine.AURORA_MYSQL, ClusterEngine.AURORA_POSTGRESQL]


def get_object_value(obj: Any, attr: str) -> Any:
    """Retrieves an arbitrary attribute value from an object.

    Nested attributes can be specified using dot notation,
    e.g. 'parent.child'.

    :param object obj:
        A valid Python object.
    :param str attr:
        The attribute name of the value to retrieve from the object.
    :returns:
        The attribute value, if it exists, or None.
    :rtype:
        any
    """
    keys = attr.split(".")
    val = obj
    for key in keys:
        if hasattr(val, key):
            val = getattr(val, key)
        else:
            return None
    return val


def merge_filters(
    filters_to_update: Optional[Dict[str, Any]], filters_to_merge: Dict[str, Any]
) -> Dict[str, Any]:
    """Given two groups of filters, merge the second into the first.

    List values are appended instead of overwritten:

    >>> merge_filters({'filter-name': ['value1']}, {'filter-name':['value2']})
    >>> {'filter-name': ['value1', 'value2']}

    :param filters_to_update:
        The filters to update.
    :type filters_to_update:
        dict[str, list] or None
    :param filters_to_merge:
        The filters to merge.
    :type filters_to_merge:
        dict[str, list] or None
    :returns:
        The updated filters.
    :rtype:
        dict[str, list]
    """
    if filters_to_update is None:
        filters_to_update = {}
    if filters_to_merge is None:
        filters_to_merge = {}
    merge_dicts(filters_to_update, filters_to_merge, append_lists=True)
    return filters_to_update


def validate_filters(
    filters: Dict[str, Any], filter_defs: Dict[str, FilterDef]
) -> None:
    """Validates filters against a set of filter definitions.

    Raises standard Python exceptions which should be caught
    and translated to an appropriate AWS/Moto exception higher
    up the call stack.

    :param dict[str, list] filters:
        The filters to validate.
    :param dict[str, FilterDef] filter_defs:
        The filter definitions to validate against.
    :returns: None
    :rtype: None
    :raises KeyError:
        if filter name not found in the filter definitions.
    :raises ValueError:
        if filter values is an empty list.
    :raises NotImplementedError:
        if `moto` does not yet support this filter.
    """
    for filter_name, filter_values in filters.items():
        filter_def = filter_defs.get(filter_name)
        if filter_def is None:
            raise KeyError(f"Unrecognized filter name: {filter_name}")
        if not filter_values:
            raise ValueError(f"The list of {filter_def.description} must not be empty.")
        if filter_def.attrs_to_check is None:
            raise NotImplementedError(
                f"{filter_name} filter has not been implemented in Moto yet."
            )


def apply_filter(resources: Any, filters: Any, filter_defs: Any) -> Any:
    """Apply an arbitrary filter to a group of resources.

    :param dict[str, object] resources:
        A dictionary mapping resource identifiers to resource objects.
    :param dict[str, list] filters:
        The filters to apply.
    :param dict[str, FilterDef] filter_defs:
        The supported filter definitions for the resource type.
    :returns:
        The filtered collection of resources.
    :rtype:
        dict[str, object]
    """
    resources_filtered = OrderedDict()
    for identifier, obj in resources.items():
        matches_filter = False
        for filter_name, filter_values in filters.items():
            filter_def = filter_defs.get(filter_name)
            for attr in filter_def.attrs_to_check:
                if get_object_value(obj, attr) in filter_values:
                    matches_filter = True
                    break
            else:
                matches_filter = False
            if not matches_filter:
                break
        if matches_filter:
            resources_filtered[identifier] = obj
    return resources_filtered


def get_start_date_end_date(
    base_date: str, window: str
) -> Tuple[datetime.datetime, datetime.datetime]:
    """Gets the start date and end date given DDD:HH24:MM-DDD:HH24:MM.

    :param base_date:
        type datetime
    :param window:
        DDD:HH24:MM-DDD:HH24:MM
    :returns:
        Start and End Date in datetime format
    :rtype:
        tuple
    """
    days = {"mon": 1, "tue": 2, "wed": 3, "thu": 4, "fri": 5, "sat": 6, "sun": 7}
    start = datetime.datetime.strptime(
        base_date + " " + window[4:9], "%d/%m/%y %H:%M"
    ) + datetime.timedelta(days=days[window[0:3]])
    end = datetime.datetime.strptime(
        base_date + " " + window[14::], "%d/%m/%y %H:%M"
    ) + datetime.timedelta(days=days[window[10:13]])
    return start, end


def get_start_date_end_date_from_time(
    base_date: str, window: str
) -> Tuple[datetime.datetime, datetime.datetime, bool]:
    """Gets the start date and end date given HH24:MM-HH24:MM.

    :param window:
        HH24:MM-HH24:MM
    :returns:
        Start and End Date in datetime format
        along with flag for spills over a day
        This is useful when determine time overlaps
    :rtype:
        tuple
    """
    times = window.split("-")
    spillover = False
    start = datetime.datetime.strptime(base_date + " " + times[0], "%d/%m/%y %H:%M")
    end = datetime.datetime.strptime(base_date + " " + times[1], "%d/%m/%y %H:%M")
    if end < start:
        end += datetime.timedelta(days=1)
        spillover = True
    return start, end, spillover


def get_overlap_between_two_date_ranges(
    start_time_1: datetime.datetime,
    end_time_1: datetime.datetime,
    start_time_2: datetime.datetime,
    end_time_2: datetime.datetime,
) -> int:
    """
    Determines overlap between 2 date ranges. Returns the overlap in seconds.
    """
    latest_start = max(start_time_1, start_time_2)
    earliest_end = min(end_time_1, end_time_2)
    delta = earliest_end - latest_start
    return (delta.days * SECONDS_IN_ONE_DAY) + delta.seconds


def valid_preferred_maintenance_window(
    maintenance_window: Any, backup_window: Any
) -> Optional[str]:
    """Determines validity of preferred_maintenance_window

    :param maintenance_windown:
        type DDD:HH24:MM-DDD:HH24:MM
    :param backup_window:
        type HH24:MM-HH24:MM
    :returns:
        message
    :rtype:
        str
    """
    MINUTES_30 = 1800
    HOURS_24 = 86400
    base_date = datetime.datetime.now().strftime("%d/%m/%y")
    try:
        p = re.compile(
            "([a-z]{3}):([0-9]{2}):([0-9]{2})-([a-z]{3}):([0-9]{2}):([0-9]{2})"
        )
        if len(maintenance_window) != 19 or re.search(p, maintenance_window) is None:
            return f"Invalid maintenance window format: {maintenance_window}. Should be specified as a range ddd:hh24:mi-ddd:hh24:mi (24H Clock UTC). Example: Sun:23:45-Mon:00:15"
        if backup_window:
            (
                backup_window_start,
                backup_window_end,
                backup_spill,
            ) = get_start_date_end_date_from_time(base_date, backup_window)
            (
                maintenance_window_start,
                maintenance_window_end,
                maintenance_spill,
            ) = get_start_date_end_date_from_time(
                base_date, maintenance_window[4:10] + maintenance_window[14::]
            )
            if (
                get_overlap_between_two_date_ranges(
                    backup_window_start,
                    backup_window_end,
                    maintenance_window_start,
                    maintenance_window_end,
                )
                >= 0
            ):
                return "The backup window and maintenance window must not overlap."

            # Due to spill overs, adjust the windows
            elif maintenance_spill:
                backup_window_start += datetime.timedelta(days=1)
                backup_window_end += datetime.timedelta(days=1)
            elif backup_spill:
                maintenance_window_start += datetime.timedelta(days=1)
                maintenance_window_end += datetime.timedelta(days=1)

            # If spills, rerun overlap test with adjusted windows
            if maintenance_spill or backup_spill:
                if (
                    get_overlap_between_two_date_ranges(
                        backup_window_start,
                        backup_window_end,
                        maintenance_window_start,
                        maintenance_window_end,
                    )
                    >= 0
                ):
                    return "The backup window and maintenance window must not overlap."

        maintenance_window_start, maintenance_window_end = get_start_date_end_date(
            base_date, maintenance_window
        )
        delta = maintenance_window_end - maintenance_window_start
        delta_seconds = delta.seconds + (delta.days * SECONDS_IN_ONE_DAY)
        if delta_seconds >= MINUTES_30 and delta_seconds <= HOURS_24:
            return None
        elif delta_seconds >= 0 and delta_seconds <= MINUTES_30:
            return "The maintenance window must be at least 30 minutes."
        else:
            return "Maintenance window must be less than 24 hours."
    except Exception:
        return f"Invalid day:hour:minute value: {maintenance_window}"


ORDERABLE_DB_INSTANCE_ENCODING = {
    "Engine": "E",
    "EngineVersion": "EV",
    "DBInstanceClass": "DBIC",
    "LicenseModel": "L",
    "AvailabilityZones": "AZ",
    "MultiAZCapable": "MC",
    "ReadReplicaCapable": "RC",
    "Vpc": "V",
    "SupportsStorageEncryption": "SE",
    "StorageType": "ST",
    "SupportsIops": "SI",
    "SupportsEnhancedMonitoring": "SM",
    "SupportsIAMDatabaseAuthentication": "SIAM",
    "SupportsPerformanceInsights": "SPI",
    "AvailableProcessorFeatures": "APF",
    "SupportedEngineModes": "SEM",
    "SupportsKerberosAuthentication": "SK",
    "OutpostCapable": "O",
    "SupportedActivityStreamModes": "SSM",
    "SupportsGlobalDatabases": "SGD",
    "SupportsClusters": "SC",
    "SupportedNetworkTypes": "SN",
    "SupportsStorageThroughput": "SST",
}
ORDERABLE_DB_INSTANCE_DECODING = {
    v: k for (k, v) in ORDERABLE_DB_INSTANCE_ENCODING.items()
}


def encode_orderable_db_instance(db: Dict[str, Any]) -> Dict[str, Any]:
    encoded = copy.deepcopy(db)
    if "AvailabilityZones" in encoded:
        encoded["AvailabilityZones"] = [
            az["Name"] for az in encoded["AvailabilityZones"]
        ]
    return {
        ORDERABLE_DB_INSTANCE_ENCODING.get(key, key): value
        for key, value in encoded.items()
    }


def decode_orderable_db_instance(db: Dict[str, Any]) -> Dict[str, Any]:
    decoded = copy.deepcopy(db)
    decoded_az = ORDERABLE_DB_INSTANCE_ENCODING.get(
        "AvailabilityZones", "AvailabilityZones"
    )
    if decoded_az in decoded:
        decoded["AvailabilityZones"] = [{"Name": az} for az in decoded[decoded_az]]
    return {
        ORDERABLE_DB_INSTANCE_DECODING.get(key, key): value
        for key, value in decoded.items()
    }

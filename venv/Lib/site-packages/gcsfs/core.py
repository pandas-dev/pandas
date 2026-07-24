"""
Google Cloud Storage pythonic interface
"""

import asyncio
import io
import json
import logging
import mimetypes
import os
import posixpath
import re
import sys
import uuid
import warnings
import weakref
from datetime import datetime, timedelta
from glob import has_magic
from urllib.parse import parse_qs
from urllib.parse import quote as quote_urllib
from urllib.parse import urlsplit

import aiohttp
import fsspec
from fsspec import asyn
from fsspec.callbacks import NoOpCallback
from fsspec.implementations.http import get_client
from fsspec.utils import other_paths, setup_logging, stringify_path

from . import __version__ as version
from ._dircache import DirCacheUpdater
from .checkers import get_consistency_checker
from .concurrency import parallel_tasks_first_completed, split_range
from .credentials import GoogleCredentials
from .inventory_report import InventoryReport
from .retry import errs, retry_request, validate_response
from .zb_hns_utils import DEFAULT_CONCURRENCY, MAX_PREFETCH_SIZE

logger = logging.getLogger("gcsfs")


if "GCSFS_DEBUG" in os.environ:
    setup_logging(logger=logger, level=os.getenv("GCSFS_DEBUG"))


# client created 2018-01-16
ACLs = {
    "authenticatedread",
    "bucketownerfullcontrol",
    "bucketownerread",
    "private",
    "projectprivate",
    "publicread",
}
bACLs = {
    "authenticatedRead",
    "private",
    "projectPrivate",
    "publicRead",
    "publicReadWrite",
}
DEFAULT_PROJECT = os.getenv("GCSFS_DEFAULT_PROJECT", "")

GCS_MIN_BLOCK_SIZE = 2**18
GCS_MAX_BLOCK_SIZE = 2**28
DEFAULT_BLOCK_SIZE = 5 * 2**20

SUPPORTED_FIXED_KEY_METADATA = {
    "content_encoding": "contentEncoding",
    "cache_control": "cacheControl",
    "content_disposition": "contentDisposition",
    "content_language": "contentLanguage",
    "custom_time": "customTime",
}

# Define allowed parameters for the GCS list API
_VALID_LIST_PARAMS = {
    "delimiter",
    "prefix",
    "startOffset",
    "endOffset",
    "maxResults",
    "versions",
    "pageToken",
    "includeFoldersAsPrefixes",
}


def quote(s):
    """
    Quote characters to be safe for URL paths.
    Also quotes '/'.

    Parameters
    ----------
    s: input URL/portion

    Returns
    -------
    corrected URL
    """
    # Encode everything, including slashes
    return quote_urllib(s, safe="")


def norm_path(path):
    """
    Canonicalize path to '{bucket}/{name}' form.

    Used by petastorm, do not remove.
    """
    bucket, name, _ = GCSFileSystem._split_path(path)
    return "/".join((bucket, name))


async def _req_to_text(r):
    async with r:
        return (await r.read()).decode()


class UnclosableBytesIO(io.BytesIO):
    """Prevent closing BytesIO to avoid errors during retries."""

    def close(self):
        """Reset stream position for next retry."""
        self.seek(0)


def _gcp_universe_domain():
    return os.getenv("GOOGLE_CLOUD_UNIVERSE_DOMAIN", "googleapis.com")


def _location():
    """
    Resolves GCS HTTP location as http[s]://host

    Enables storage emulation for integration tests.

    Returns
    -------
    valid http location
    """
    _emulator_location = os.getenv("STORAGE_EMULATOR_HOST", "")
    if _emulator_location not in {"default", "", None}:
        if not any(
            _emulator_location.startswith(scheme) for scheme in ("http://", "https://")
        ):
            _emulator_location = f"http://{_emulator_location}"
        return _emulator_location

    return f"https://storage.{_gcp_universe_domain()}"


def _chunks(lst, n):
    """
    Yield evenly-sized chunks from a list.

    Implementation based on https://stackoverflow.com/a/312464.
    """
    for i in range(0, len(lst), n):
        yield lst[i : i + n]


def _coalesce_generation(*args):
    """Helper to coalesce a list of object generations down to one."""
    generations = set(args)
    if None in generations:
        generations.remove(None)
    if len(generations) > 1:
        raise ValueError(
            "Cannot coalesce generations where more than one are defined,"
            f" {generations}"
        )
    elif len(generations) == 0:
        return None
    else:
        return generations.pop()


def _is_directory_marker(entry):
    return entry["size"] == 0 and entry["name"].endswith("/")


class GCSFileSystem(DirCacheUpdater, asyn.AsyncFileSystem):
    r"""
    Connect to Google Cloud Storage.

    The following modes of authentication are supported:

    - ``token=None``, GCSFS will attempt to guess your credentials in the
      following order: gcloud CLI default, gcsfs cached token, google compute
      metadata service, anonymous.
    - ``token='google_default'``, your default gcloud credentials will be used,
      which are typically established by doing ``gcloud login`` in a terminal.
    - ``token='cache'``, credentials from previously successful gcsfs
      authentication will be used (use this after "browser" auth succeeded)
    - ``token='anon'``, no authentication is performed, and you can only
      access data which is accessible to allUsers (in this case, the project and
      access level parameters are meaningless)
    - ``token='browser'``, you get an access code with which you can
      authenticate via a specially provided URL
    - if ``token='cloud'``, we assume we are running within google compute
      or google container engine, and query the internal metadata directly for
      a token.
    - you may supply a token generated by the
      [gcloud](https://cloud.google.com/sdk/docs/)
      utility; this is either a python dictionary, the name of a file
      containing the JSON returned by logging in with the gcloud CLI tool,
      or a Credentials object. gcloud typically stores its tokens in locations
      such as
      ``~/.config/gcloud/application_default_credentials.json``,
      ``~/.config/gcloud/credentials``, or
      ``~\AppData\Roaming\gcloud\credentials``, etc.

    Specific methods, (eg. ``ls``, ``info``, ...) may return object details from GCS.
    These detailed listings include the
    [object resource](https://cloud.google.com/storage/docs/json_api/v1/objects#resource)

    GCS *does not* include  "directory" objects but instead generates
    directories by splitting
    [object names](https://cloud.google.com/storage/docs/key-terms).
    This means that, for example,
    a directory does not need to exist for an object to be created within it.
    Creating an object implicitly creates it's parent directories, and removing
    all objects from a directory implicitly deletes the empty directory.

    `GCSFileSystem` generates listing entries for these implied directories in
    listing apis with the  object properties:

        - "name" : string
            The "{bucket}/{name}" path of the dir, used in calls to
            GCSFileSystem or GCSFile.
        - "bucket" : string
            The name of the bucket containing this object.
        - "kind" : 'storage#object'
        - "size" : 0
        - "storageClass" : 'DIRECTORY'
        - type: 'directory' (fsspec compat)

    GCSFileSystem maintains a per-implied-directory cache of object listings and
    fulfills all object information and listing requests from cache. This implied, for example, that objects
    created via other processes *will not* be visible to the GCSFileSystem until the cache
    refreshed. Calls to GCSFileSystem.open and calls to GCSFile are not affected by this cache.

    Note that directory listings are cached by default, because fetching those listings can be expensive. This is
    contrary to local filesystem behaviour. The cache will be cleared if writing from this instance, but it can
    become stale and return incorrect results if the storage is written to from another process/machine.
    If you anticipate this possibility, you can set the use_listings_cache and listings_expiry_time arguments
    to configure the caching, call `.invalidate_cache()` when required, or pass `refresh=True` to the
    various listing methods.

    In the default case the cache is never expired. This may be controlled via the ``cache_timeout``
    GCSFileSystem parameter or via explicit calls to ``GCSFileSystem.invalidate_cache``.

    NOTE on "exclusive" mode: mode=="create"" (in pipe and put) and open(mode="xb") are supported on an
    experimental basis. The test harness does not currently support this, so use at your
    own risk.

    Parameters
    ----------
    project : string
        project_id to work under. Note that this is not the same as, but often
        very similar to, the project name.
        This is required in order
        to list all the buckets you have access to within a project and to
        create/delete buckets, or update their access policies.
        If ``token='google_default'``, the value is overridden by the default,
        if ``token='anon'``, the value is ignored.
    access : one of {'read_only', 'read_write', 'full_control'}
        Full control implies read/write as well as modifying metadata,
        e.g., access control.
    token: None, dict or string
        (see description of authentication methods, above)
    consistency: 'none', 'size', 'md5'
        Check method when writing files. Can be overridden in open().
    cache_timeout: float, seconds
        Cache expiration time in seconds for object metadata cache.
        Set cache_timeout <= 0 for no caching, None for no cache expiration.
    secure_serialize: bool (deprecated)
    requester_pays : bool, or str default False
        Whether to use requester-pays requests. This will include your
        project ID `project` in requests as the `userProject`, and you'll be
        billed for accessing data from requester-pays buckets. Optionally,
        pass a project-id here as a string to use that as the `userProject`.
    session_kwargs: dict
        passed on to ``aiohttp.ClientSession``; can contain, for example,
        proxy settings.
    endpoint_url: str
        If given, use this URL (format protocol://host:port , *without* any
        path part) for communication. If not given, defaults to the value
        of environment variable "STORAGE_EMULATOR_HOST"; if that is not set
        either, will use the standard Google endpoint.
    default_location: str
        Default location where buckets are created, like 'US' or 'EUROPE-WEST3'.
        You can find a list of all available locations here:
        https://cloud.google.com/storage/docs/locations#available-locations
    version_aware: bool
        Whether to support object versioning. If enabled this will require the
        user to have the necessary permissions for dealing with versioned objects.
    """

    scopes = {"read_only", "read_write", "full_control"}
    retries = 6  # number of retries on http failure
    default_block_size = DEFAULT_BLOCK_SIZE
    protocol = "gs", "gcs"
    async_impl = True
    MIN_CHUNK_SIZE_FOR_CONCURRENCY = 5 * 1024 * 1024

    def __init__(
        self,
        project=DEFAULT_PROJECT,
        access="full_control",
        token=None,
        block_size=None,
        consistency="none",
        cache_timeout=None,
        secure_serialize=True,
        check_connection=None,
        requests_timeout=None,
        requester_pays=False,
        asynchronous=False,
        session_kwargs=None,
        loop=None,
        timeout=None,
        endpoint_url=None,
        default_location=None,
        version_aware=False,
        **kwargs,
    ):
        if cache_timeout is not None:
            kwargs["listings_expiry_time"] = cache_timeout
        super().__init__(
            self,
            asynchronous=asynchronous,
            loop=loop,
            **kwargs,
        )
        if access not in self.scopes:
            raise ValueError(f"access must be one of {self.scopes}")
        if project is None:
            warnings.warn("GCS project not set - cannot list or create buckets")
        if block_size is not None:
            self.default_block_size = block_size
        self.requester_pays = requester_pays
        self.consistency = consistency
        self.cache_timeout = cache_timeout or kwargs.pop("listings_expiry_time", None)
        self.requests_timeout = requests_timeout
        self.timeout = timeout
        self._session = None
        self._endpoint = endpoint_url
        self.session_kwargs = session_kwargs or {}
        self.default_location = default_location
        self.version_aware = version_aware

        if check_connection:
            warnings.warn(
                "The `check_connection` argument is deprecated and will be removed in a future release.",
                DeprecationWarning,
            )

        self.credentials = GoogleCredentials(
            project, access, token, on_google=self.on_google
        )

    @property
    def _location(self):
        return self._endpoint or _location()

    @property
    def base(self):
        return f"{self._location}/storage/v1/"

    @property
    def batch_url_base(self):
        return f"{self._location}/batch/storage/v1"

    @property
    def project(self):
        return self.credentials.project

    # This threshold applies to the standard bucket, whereas the zonal bucket
    # uses a 5MB threshold. This difference exists because the standard bucket
    # lacks the `DirectMemmoveBuffer` implementation used in the zonal bucket.
    async def _get_threshold_for_disk_reads(self, bucket):
        return 100 * 1024 * 1024

    # Clean up the aiohttp session
    #
    # This can run from the main thread if invoked via the weakref callback.
    # This can happen even if the `loop` parameter belongs to another thread
    # (e.g. the fsspec IO worker). The control flow here is intended to attempt
    # in-thread asynchronous cleanup first, then fallback to synchronous
    # cleanup (which can handle cross-thread calls).
    @staticmethod
    def close_session(loop, session: aiohttp.ClientSession, asynchronous=False):
        if session is None or session.closed:
            return
        force_close = False
        try:
            current_loop = asyncio.get_running_loop()
        except RuntimeError:
            current_loop = None
        if loop:
            # an explicit loop was set
            if loop.is_running():
                loop.create_task(session.close())
            else:
                force_close = True
        elif current_loop is not None and current_loop.is_running() and asynchronous:
            # running in a concurrnet context
            current_loop.create_task(session.close())
        elif asyn.loop[0] is not None and asyn.loop[0].is_running():
            try:
                asyn.sync(asyn.loop[0], session.close, timeout=0.1)
            except fsspec.FSTimeoutError:
                force_close = True
        else:
            force_close = True
        if force_close:
            # during shutdown, this is the fallback
            connector = getattr(session, "_connector", None)
            if connector is not None:
                # close after loop is dead
                connector._close()

    async def _set_session(self):
        if self._session is None:
            self._session = await get_client(**self.session_kwargs)
            weakref.finalize(
                self, self.close_session, self.loop, self._session, self.asynchronous
            )
        return self._session

    @property
    def session(self):
        if self.asynchronous and self._session is None:
            raise RuntimeError("Please await _connect* before anything else")
        return self._session

    @classmethod
    def _strip_protocol(cls, path):
        if isinstance(path, list):
            return [cls._strip_protocol(p) for p in path]
        path = stringify_path(path)
        protos = (cls.protocol,) if isinstance(cls.protocol, str) else cls.protocol
        for protocol in protos:
            if path.startswith(protocol + "://"):
                path = path[len(protocol) + 3 :]
            elif path.startswith(protocol + "::"):
                path = path[len(protocol) + 2 :]
        # use of root_marker to make minimum required path, e.g., "/"
        return path or cls.root_marker

    @classmethod
    def _get_kwargs_from_urls(cls, path):
        _, _, generation = cls._split_path(path, version_aware=True)
        if generation is not None:
            return {"version_aware": True}
        return {}

    def _get_params(self, kwargs):
        params = {k: v for k, v in kwargs.items() if v is not None}
        # needed for requester pays buckets
        if self.requester_pays:
            if isinstance(self.requester_pays, str):
                user_project = self.requester_pays
            else:
                user_project = self.project
            params["userProject"] = user_project
        return params

    def _get_headers(self, headers, cache_type=None):
        out = {}
        if headers is not None:
            out.update(headers)
        if "User-Agent" not in out:
            ua = "python-gcsfs/" + version
            if cache_type:
                ua += f" cache_type/{cache_type}"
            out["User-Agent"] = ua
        self.credentials.apply(out)
        return out

    def _format_path(self, path, args):
        if not path.startswith("http"):
            path = self.base + path

        if args:
            path = path.format(*[quote(p) for p in args])
        return path

    @retry_request(retries=retries)
    async def _request(
        self,
        method,
        path,
        *args,
        headers=None,
        json=None,
        data=None,
        cache_type=None,
        **kwargs,
    ):
        await self._set_session()
        if hasattr(data, "seek"):
            data.seek(0)
        async with self.session.request(
            method=method,
            url=self._format_path(path, args),
            params=self._get_params(kwargs),
            json=json,
            headers=self._get_headers(headers, cache_type=cache_type),
            data=data,
            timeout=self.requests_timeout,
        ) as r:
            status = r.status
            headers = r.headers
            info = r.request_info  # for debug only
            contents = await r.read()

            validate_response(status, contents, path, args)
            return status, headers, info, contents

    async def _call(
        self, method, path, *args, json_out=False, info_out=False, **kwargs
    ):
        logger.debug(f"{method.upper()}: {path}, {args}, {kwargs.get('headers')}")
        status, headers, info, contents = await self._request(
            method, path, *args, **kwargs
        )
        if json_out:
            return json.loads(contents)
        elif info_out:
            return info
        else:
            return headers, contents

    call = asyn.sync_wrapper(_call)

    @property
    def buckets(self):
        """Return list of available project buckets."""
        return [
            b["name"]
            for b in asyn.sync(self.loop, self._list_buckets, timeout=self.timeout)
        ]

    def _process_object(self, bucket, object_metadata):
        """Process object resource into gcsfs object information format.

        Process GCS object resource via type casting and attribute updates to
        the cache-able gcsfs object information format. Returns an updated copy
        of the object resource.

        (See https://cloud.google.com/storage/docs/json_api/v1/objects#resource)
        """
        result = dict(object_metadata)
        result["size"] = int(object_metadata.get("size", 0))
        result["name"] = posixpath.join(bucket, object_metadata["name"])
        result["type"] = "file"
        # Translate time metadata from GCS names to fsspec standard names.
        # TODO(issues/559): Remove legacy names `updated` and `timeCreated`?
        if "updated" in object_metadata:
            result["mtime"] = self._parse_timestamp(object_metadata["updated"])
        if "timeCreated" in object_metadata:
            result["ctime"] = self._parse_timestamp(object_metadata["timeCreated"])
        if "generation" in object_metadata or "metageneration" in object_metadata:
            result["generation"] = object_metadata.get("generation")
            result["metageneration"] = object_metadata.get("metageneration")

        return result

    async def _make_bucket_requester_pays(self, path, state=True):
        # this is really some form of setACL/chmod
        # perhaps should be automatic if gcs.requester_pays
        json = {"billing": {"requesterPays": state}}
        await self._call("PATCH", f"b/{path}", json=json)

    make_bucket_requester_pays = asyn.sync_wrapper(_make_bucket_requester_pays)

    async def _get_object(self, path):
        """Return object information at the given path."""
        bucket, key, generation = self.split_path(path)

        # Check if parent dir is in listing cache
        listing = self._ls_from_cache(path)
        if listing:
            name = "/".join((bucket, key))
            for file_details in listing:
                if (
                    file_details["type"] == "file"
                    and file_details["name"] == name
                    and (not generation or file_details.get("generation") == generation)
                ):
                    return file_details
            else:
                raise FileNotFoundError(path)

        if not key:
            # Attempt to "get" the bucket root, return error instead of
            # listing.
            raise FileNotFoundError(path)

        res = None
        # Work around various permission settings. Prefer an object get (storage.objects.get), but
        # fall back to a bucket list + filter to object name (storage.objects.list).
        try:
            res = await self._call(
                "GET", "b/{}/o/{}", bucket, key, json_out=True, generation=generation
            )
        except OSError as e:
            if not str(e).startswith("Forbidden"):
                raise
            resp = await self._call(
                "GET",
                "b/{}/o",
                bucket,
                json_out=True,
                prefix=key,
                maxResults=1 if not generation else None,
                versions="true" if generation else None,
            )
            for item in resp.get("items", []):
                if item["name"] == key and (
                    not generation or item.get("generation") == generation
                ):
                    res = item
                    break
            if res is None:
                raise FileNotFoundError(path)
        return self._process_object(bucket, res)

    async def _list_objects(self, path, prefix="", versions=False, **kwargs):
        bucket, key, generation = self.split_path(path)
        path = path.rstrip("/")

        # NOTE: the inventory report logic is experimental.
        inventory_report_info = kwargs.get("inventory_report_info", None)

        # Only attempt to list from the cache when the user does not use
        # the inventory report service.
        if not inventory_report_info:
            try:
                clisting = self._ls_from_cache(path)
                hassubdirs = clisting and any(
                    c["name"].rstrip("/") == path and c["type"] == "directory"
                    for c in clisting
                )
                if clisting and not hassubdirs:
                    return clisting
            except FileNotFoundError:
                # not finding a bucket in list of "my" buckets is OK
                if key:
                    raise

        items, prefixes = await self._do_list_objects(
            path,
            prefix=prefix,
            versions=versions,
            **kwargs,
        )

        pseudodirs = [
            {
                "bucket": bucket,
                "name": bucket + "/" + prefix.strip("/"),
                "size": 0,
                "storageClass": "DIRECTORY",
                "type": "directory",
            }
            for prefix in prefixes
        ]
        if not (items + pseudodirs):
            if key:
                return [await self._get_object(path)]
            else:
                return []
        out = pseudodirs + items

        use_snapshot_listing = inventory_report_info and inventory_report_info.get(
            "use_snapshot_listing"
        )

        max_results = kwargs.get("max_results")

        # Don't cache prefixed/partial listings, in addition to
        # not using the inventory report service to do listing directly.
        if not prefix and not use_snapshot_listing and not max_results:
            self.dircache[path] = out
        return out

    async def _do_list_objects(
        self,
        path,
        max_results=None,
        delimiter="/",
        prefix="",
        versions=False,
        **kwargs,
    ):
        """Object listing for the given {bucket}/{prefix}/ path."""
        bucket, _path, generation = self.split_path(path)
        _path = "" if not _path else _path.rstrip("/") + "/"
        prefix = f"{_path}{prefix}" or None
        versions = bool(versions or generation)

        # Page size of 5000 is officially supported across GCS.
        default_page_size = 5000

        # NOTE: the inventory report logic is experimental.
        inventory_report_info = kwargs.get("inventory_report_info", None)

        # Check if the user has configured inventory report option.
        if inventory_report_info is not None:
            items, prefixes = await InventoryReport.fetch_snapshot(
                gcs_file_system=self,
                inventory_report_info=inventory_report_info,
                prefix=prefix,
            )

            use_snapshot_listing = inventory_report_info.get("use_snapshot_listing")

            # If the user wants to rely on the snapshot from the inventory report
            # for listing, directly return the results.
            if use_snapshot_listing:
                return items, prefixes

            # Otherwise, use the snapshot to initiate concurrent listing.
            return await self._concurrent_list_objects_helper(
                items=items,
                bucket=bucket,
                delimiter=delimiter,
                prefix=prefix,
                versions=versions,
                page_size=default_page_size,
                **kwargs,
            )

        # If the user has not configured inventory report, proceed to use
        # sequential listing.
        else:
            return await self._sequential_list_objects_helper(
                bucket=bucket,
                delimiter=delimiter,
                start_offset=None,
                end_offset=None,
                prefix=prefix,
                versions=versions,
                max_results=max_results,
                **kwargs,
            )

    async def _concurrent_list_objects_helper(
        self, items, bucket, delimiter, prefix, versions, page_size, **kwargs
    ):
        """
        Lists objects using coroutines, using the object names from the inventory
        report to split up the ranges.
        """

        # Extract out the names of the objects fetched from the inventory report.
        snapshot_object_names = sorted([item["name"] for item in items])

        # Determine the number of coroutines needed to concurrent listing.
        # Ideally, want each coroutine to fetch a single page of objects.
        num_coroutines = len(snapshot_object_names) // page_size + 1
        num_objects_per_coroutine = len(snapshot_object_names) // num_coroutines

        start_offsets = []
        end_offsets = []

        # Calculate the split splits of each coroutine (start offset and end offset).
        for i in range(num_coroutines):
            range_start = i * num_objects_per_coroutine
            if i == num_coroutines - 1:
                range_end = len(snapshot_object_names)
            else:
                range_end = range_start + num_objects_per_coroutine

            if range_start == 0:
                prefix_start = None
            else:
                prefix_start = snapshot_object_names[range_start]

            if range_end == len(snapshot_object_names):
                prefix_end = None
            else:
                prefix_end = snapshot_object_names[range_end]

            start_offsets.append(prefix_start)
            end_offsets.append(prefix_end)

        # Assign the coroutine all at once, and wait for them to finish listing.
        results = await asyncio.gather(
            *[
                self._sequential_list_objects_helper(
                    bucket=bucket,
                    delimiter=delimiter,
                    start_offset=start_offsets[i],
                    end_offset=end_offsets[i],
                    prefix=prefix,
                    versions=versions,
                    max_results=page_size,
                    **kwargs,
                )
                for i in range(0, len(start_offsets))
            ]
        )

        items = []
        prefixes = []

        # Concatenate the items and prefixes from each coroutine for final results.
        for i in range(len(results)):
            items_from_process, prefixes_from_process = results[i]
            items.extend(items_from_process)
            prefixes.extend(prefixes_from_process)

        return items, prefixes

    async def _sequential_list_objects_helper(
        self,
        bucket,
        delimiter,
        start_offset,
        end_offset,
        prefix,
        versions,
        max_results,
        items_per_call=1000,
        **kwargs,
    ):
        """
        Sequential list objects within the start and end offset range.
        """
        max_results = max_results if max_results else 10_000_000
        prefixes = []
        items = []
        num_items = min(items_per_call, max_results, 1000)
        page = await self._call_list_objects(
            bucket,
            delimiter=delimiter,
            prefix=prefix,
            startOffset=start_offset,
            endOffset=end_offset,
            maxResults=num_items,
            versions="true" if versions else None,
            **kwargs,
        )

        prefixes.extend(page.get("prefixes", []))
        items.extend(page.get("items", []))
        next_page_token = page.get("nextPageToken", None)

        while len(items) + len(prefixes) < max_results and next_page_token is not None:
            num_items = min(
                items_per_call, max_results - (len(items) + len(prefixes)), 1000
            )
            page = await self._call_list_objects(
                bucket,
                delimiter=delimiter,
                prefix=prefix,
                startOffset=start_offset,
                endOffset=end_offset,
                maxResults=num_items,
                pageToken=next_page_token,
                versions="true" if versions else None,
                **kwargs,
            )

            assert page["kind"] == "storage#objects"
            prefixes.extend(page.get("prefixes", []))
            items.extend(page.get("items", []))
            next_page_token = page.get("nextPageToken", None)

        items = [self._process_object(bucket, i) for i in items]

        return items, prefixes

    async def _call_list_objects(self, bucket, **kwargs):
        """
        Helper method to fetch a single page of object listing.
        Extracts valid GCS parameters from kwargs to prevent parameter pollution.
        """

        # Only pass valid parameters to the API call
        valid_kwargs = {k: v for k, v in kwargs.items() if k in _VALID_LIST_PARAMS}

        return await self._call(
            "GET",
            "b/{}/o",
            bucket,
            json_out=True,
            **valid_kwargs,
        )

    async def _list_buckets(self):
        """Return list of all buckets under the current project."""
        if "" not in self.dircache:
            items = []
            page = await self._call("GET", "b", project=self.project, json_out=True)

            assert page["kind"] == "storage#buckets"
            items.extend(page.get("items", []))
            next_page_token = page.get("nextPageToken", None)

            while next_page_token is not None:
                page = await self._call(
                    "GET",
                    "b",
                    project=self.project,
                    pageToken=next_page_token,
                    json_out=True,
                )

                assert page["kind"] == "storage#buckets"
                items.extend(page.get("items", []))
                next_page_token = page.get("nextPageToken", None)

            buckets = [
                {**i, "name": i["name"] + "/", "size": 0, "type": "directory"}
                for i in items
            ]
            self.dircache[""] = buckets
            return buckets
        return self.dircache[""]

    def invalidate_cache(self, path=None):
        """
        Invalidate listing cache for given path, it is reloaded on next use.

        Parameters
        ----------
        path: string or None
            If None, clear all listings cached else listings at or under given
            path.
        """
        if path is None:
            logger.debug("invalidate_cache clearing cache")
            self.dircache.clear()
        else:
            path = self._strip_protocol(path).rstrip("/")
            if not path:
                self.dircache.pop("", None)
            while path:
                self.dircache.pop(path, None)
                path = self._parent(path)

    async def _mkdir(
        self,
        path,
        acl="projectPrivate",
        default_acl="bucketOwnerFullControl",
        location=None,
        create_parents=False,
        enable_versioning=False,
        enable_object_retention=False,
        iam_configuration=None,
        **kwargs,
    ):
        """
        New bucket

        If path is more than just a bucket, will create bucket if create_parents=True;
        otherwise is a noop. If create_parents is False and bucket does not exist,
        will produce FileNotFoundError.

        Parameters
        ----------
        path: str
            bucket name. If contains '/' (i.e., looks like subdir), will
            have no effect because GCS doesn't have real directories.
        acl: string, one of bACLs
            access for the bucket itself. See:
            https://cloud.google.com/storage/docs/access-control/lists#predefined-acl
        default_acl: str, one of ACLs
            default ACL for objects created in this bucket
        location: Optional[str]
            Location where buckets are created, like 'US' or 'EUROPE-WEST3'.
            If not provided, defaults to `self.default_location`.
            You can find a list of all available locations here:
            https://cloud.google.com/storage/docs/locations#available-locations
        create_parents: bool
            If True, creates the bucket in question, if it doesn't already exist
        enable_versioning: bool
            If True, creates the bucket in question with object versioning
            enabled.
        enable_object_retention: bool
            If True, creates the bucket in question with object retention
            permanently enabled.
        iam_configuration: dict
            If provided, sets the IAM policy for the bucket. This argument
            allows setting properties such as `{publicAccessPrevention: "enforced"}`
            and `{"uniformBucketLevelAccess": {"enabled": True}}`. If passed, `acl`
            and `default_acl` are explicitly ignored.
        **kwargs
            Additional parameters passed to the API call request body. See:
            https://cloud.google.com/storage/docs/json_api/v1/buckets/insert#request-body
            for all possible options. Pass nested parameters as dictionaries, e.g.:
            `{"autoclass": {"enabled": True}}`
        """
        bucket, object, generation = self.split_path(path)
        if bucket in ["", "/"]:
            raise ValueError("Cannot create root bucket")
        if "/" in path and create_parents and await self._exists(bucket):
            # nothing to do
            return
        if "/" in path:
            if await self._exists(bucket):
                return
            if not create_parents:
                raise FileNotFoundError(bucket)

        json_data = {"name": bucket}
        location = location or self.default_location
        if location:
            json_data["location"] = location
        if enable_versioning:
            json_data["versioning"] = {"enabled": True}
        if iam_configuration:
            json_data["iamConfiguration"] = iam_configuration
            acl = None
            default_acl = None
        if kwargs:
            json_data.update(kwargs)

        await self._call(
            method="POST",
            path="b",
            predefinedAcl=acl,
            project=self.project,
            predefinedDefaultObjectAcl=default_acl,
            enableObjectRetention=str(enable_object_retention).lower(),
            json=json_data,
            json_out=True,
        )
        self.invalidate_cache(bucket)
        self.invalidate_cache("")

    mkdir = asyn.sync_wrapper(_mkdir)

    async def _rmdir(self, bucket):
        """Delete an empty bucket

        Parameters
        ----------
        bucket: str
            bucket name. If contains '/' (i.e., looks like subdir), will
            have no effect because GCS doesn't have real directories.
        """
        bucket = bucket.rstrip("/")
        if "/" in bucket:
            return
        await self._call("DELETE", "b/" + bucket, json_out=False)
        self.invalidate_cache(bucket)
        self.invalidate_cache("")

    rmdir = asyn.sync_wrapper(_rmdir)

    def modified(self, path):
        return self.info(path)["mtime"]

    def created(self, path):
        return self.info(path)["ctime"]

    def _parse_timestamp(self, timestamp):
        assert timestamp.endswith("Z")
        timestamp = timestamp[:-1]
        timestamp = timestamp + "0" * (6 - len(timestamp.rsplit(".", 1)[-1]))
        return datetime.fromisoformat(timestamp + "+00:00")

    async def _info(self, path, generation=None, **kwargs):
        """File information about this path."""
        path = self._strip_protocol(path).rstrip("/")
        if "/" not in path:
            try:
                out = await self._call("GET", f"b/{path}", json_out=True)
                out.update(size=0, type="directory")
            except OSError:
                # GET bucket failed, try ls; will have no metadata
                exists = await self._ls(path)
                if exists:
                    out = {"name": path, "size": 0, "type": "directory"}
                else:
                    raise FileNotFoundError(path)
            return out
        # Check directory cache for parent dir
        parent_path = self._parent(path)
        parent_cache = self._ls_from_cache(parent_path)
        bucket, key, path_generation = self.split_path(path)
        generation = _coalesce_generation(generation, path_generation)
        if parent_cache:
            name = "/".join((bucket, key))
            for o in parent_cache:
                if o["name"].rstrip("/") == name and (
                    not generation or o.get("generation") == generation
                ):
                    return o
        if self._ls_from_cache(path):
            # this is a directory
            return {
                "bucket": bucket,
                "name": path,
                "size": 0,
                "storageClass": "DIRECTORY",
                "type": "directory",
            }

        async with parallel_tasks_first_completed(
            [
                self._get_object(path),
                self._get_directory_info(path, bucket, key, generation),
            ]
        ) as (tasks, done, pending):
            get_object_task, get_directory_info_task = tasks

            try:
                get_object_res = await get_object_task
                if not _is_directory_marker(get_object_res):
                    return get_object_res
            except FileNotFoundError:
                pass
            return await get_directory_info_task

    async def _get_directory_info(self, path, bucket, key, generation):
        """
        Internal method to check if a path is a directory by listing objects.
        """
        out = await self._list_objects(path, max_results=1)
        exact = next((o for o in out if o["name"].rstrip("/") == path), None)
        if exact and not _is_directory_marker(exact):
            # exact hit
            return exact
        elif out:
            # other stuff - must be a directory
            return {
                "bucket": bucket,
                "name": path,
                "size": 0,
                "storageClass": "DIRECTORY",
                "type": "directory",
            }
        else:
            raise FileNotFoundError(path)

    async def _ls(
        self, path, detail=False, prefix="", versions=False, refresh=False, **kwargs
    ):
        """List objects under the given '/{bucket}/{prefix} path."""
        path = self._strip_protocol(path).rstrip("/")

        if refresh:
            self.invalidate_cache(path)
        if path in ["/", ""]:
            out = await self._list_buckets()
        else:
            out = []
            dir_names = set()
            for entry in await self._list_objects(
                path, prefix=prefix, versions=versions, **kwargs
            ):
                if _is_directory_marker(entry):
                    entry = {
                        "bucket": entry["bucket"],
                        "name": path.rstrip("/"),
                        "size": 0,
                        "storageClass": "DIRECTORY",
                        "type": "directory",
                    }

                if entry["type"] == "directory":
                    if entry["name"] in dir_names:
                        continue
                    dir_names.add(entry["name"])

                if versions and "generation" in entry:
                    entry = entry.copy()
                    entry["name"] = f"{entry['name']}#{entry['generation']}"

                out.append(entry)

        out.sort(key=lambda e: (e["name"]))

        if detail:
            return out
        else:
            return [o["name"] for o in out]

    def url(self, path, generation=None):
        """Get HTTP URL of the given path"""
        u = "{}/download/storage/v1/b/{}/o/{}?alt=media{}"
        bucket, object, path_generation = self.split_path(path)
        generation = _coalesce_generation(generation, path_generation)
        object = quote(object)
        return u.format(
            self._location,
            bucket,
            object,
            f"&generation={generation}" if generation else "",
        )

    async def _cat_file_sequential(self, path, start=None, end=None, **kwargs):
        """Simple one-shot get of file data"""
        # if start and end are both provided and valid, but start >= end, return empty bytes
        # Otherwise, _process_limits would generate an invalid HTTP range (e.g. "bytes=5-4"
        # for start=5, end=5), causing the server to return the whole file instead of nothing.
        if start is not None and end is not None and start >= end >= 0:
            return b""

        u2 = self.url(path, generation=kwargs.get("generation"))
        if start is not None or end is not None:
            head = {"Range": await self._process_limits(path, start, end)}
        else:
            head = {}

        cache_type = kwargs.get("cache_type")
        headers, out = await self._call("GET", u2, headers=head, cache_type=cache_type)
        return out

    async def _cat_file_concurrent(
        self, path, start=None, end=None, concurrency=DEFAULT_CONCURRENCY, **kwargs
    ):
        """Concurrent fetch of file data"""
        if start is None:
            start = 0
        if end is None:
            end = (await self._info(path))["size"]
        if start >= end:
            return b""

        ranges = split_range(
            end - start, concurrency, self.MIN_CHUNK_SIZE_FOR_CONCURRENCY
        )
        if len(ranges) == 1:
            return await self._cat_file_sequential(path, start=start, end=end, **kwargs)

        tasks = []

        for relative_offset, size in ranges:
            offset = start + relative_offset
            tasks.append(
                asyncio.create_task(
                    self._cat_file_sequential(
                        path, start=offset, end=offset + size, **kwargs
                    )
                )
            )

        try:
            results = await asyncio.gather(*tasks)
            return b"".join(results)
        except BaseException as e:
            for t in tasks:
                if not t.done():
                    t.cancel()
            await asyncio.gather(*tasks, return_exceptions=True)
            raise e

    async def _cat_file(
        self, path, start=None, end=None, concurrency=DEFAULT_CONCURRENCY, **kwargs
    ):
        """Simple one-shot, or concurrent get of file data"""
        if concurrency > 1:
            return await self._cat_file_concurrent(
                path, start=start, end=end, concurrency=concurrency, **kwargs
            )

        # While we could just call _cat_file_concurrent(concurrency=1), we are choosing
        # to keep it separate because concurrency code path is still in an experimental phase.
        # Once concurrency code path is stabilized, we can remove this if-else condition.
        return await self._cat_file_sequential(path, start=start, end=end, **kwargs)

    async def _getxattr(self, path, attr):
        """Get user-defined metadata attribute"""
        meta = (await self._info(path)).get("metadata", {})
        return meta[attr]

    getxattr = asyn.sync_wrapper(_getxattr)

    async def _setxattrs(
        self,
        path,
        content_type=None,
        content_encoding=None,
        fixed_key_metadata=None,
        **kwargs,
    ):
        """Set/delete/add writable metadata attributes

        Note: uses PATCH method (update), leaving unedited keys alone.
        fake-gcs-server:latest does not seem to support this.

        Parameters
        ----------
        content_type: str
            If not None, set the content-type to this value
        content_encoding: str
            This parameter is deprecated, you may use fixed_key_metadata instead.
            If not None, set the content-encoding.
            See https://cloud.google.com/storage/docs/transcoding
        fixed_key_metadata: dict
            Google metadata, in key/value pairs, supported keys:
                - cache_control
                - content_disposition
                - content_encoding
                - content_language
                - custom_time

            More info:
            https://cloud.google.com/storage/docs/metadata#mutable
        kw_args: key-value pairs like field="value" or field=None
            value must be string to add or modify, or None to delete

        Returns
        -------
        Entire metadata after update (even if only path is passed)
        """
        i_json = {"metadata": kwargs}
        if content_type is not None:
            i_json["contentType"] = content_type
        if content_encoding is not None:
            logger.warn(
                "setxattrs: content_encoding parameter is now deprecated "
                "you may use `fixed_key_metadata` instead"
            )
            i_json["contentEncoding"] = content_encoding
        i_json.update(_convert_fixed_key_metadata(fixed_key_metadata))

        bucket, key, generation = self.split_path(path)
        o_json = await self._call(
            "PATCH",
            "b/{}/o/{}",
            bucket,
            key,
            fields="metadata",
            json=i_json,
            json_out=True,
        )
        return o_json.get("metadata", {})

    setxattrs = asyn.sync_wrapper(_setxattrs)

    async def _merge(self, path, paths, acl=None):
        """Concatenate objects within a single bucket"""
        bucket, key, generation = self.split_path(path)
        source = [{"name": self.split_path(p)[1]} for p in paths]
        await self._call(
            "POST",
            "b/{}/o/{}/compose",
            bucket,
            key,
            destinationPredefinedAcl=acl,
            headers={"Content-Type": "application/json"},
            json={
                "sourceObjects": source,
                "kind": "storage#composeRequest",
                "destination": {"name": key, "bucket": bucket},
            },
        )

    merge = asyn.sync_wrapper(_merge)

    # TODO: Add async mv method in the async.py and remove from GCSFileSystem.
    async def _mv(
        self, path1, path2, recursive=False, maxdepth=None, batch_size=None, **kwargs
    ):
        if path1 == path2:
            return

        if isinstance(path1, list) and isinstance(path2, list):
            # No need to expand paths when both source and destination
            # are provided as lists
            paths1 = path1
            paths2 = path2
        else:
            source_is_str = isinstance(path1, str)
            paths1 = await self._expand_path(
                path1, maxdepth=maxdepth, recursive=recursive
            )
            if source_is_str and (not recursive or maxdepth is not None):
                # Non-recursive glob does not move directories
                paths1 = [
                    p
                    for p in paths1
                    if not (asyn.trailing_sep(p) or await self._isdir(p))
                ]
                if not paths1:
                    return

            source_is_file = len(paths1) == 1
            dest_is_dir = isinstance(path2, str) and (
                asyn.trailing_sep(path2) or await self._isdir(path2)
            )

            exists = source_is_str and (
                (has_magic(path1) and source_is_file)
                or (
                    not has_magic(path1)
                    and dest_is_dir
                    and not asyn.trailing_sep(path1)
                )
            )
            paths2 = other_paths(
                paths1,
                path2,
                exists=exists,
                flatten=not source_is_str,
            )

        batch_size = batch_size or self.batch_size
        result = await asyn._run_coros_in_chunks(
            [self._mv_file(p1, p2, **kwargs) for p1, p2 in zip(paths1, paths2)],
            batch_size=batch_size,
            return_exceptions=True,
            nofiles=True,
        )

        for res, p1 in zip(result, paths1):
            if isinstance(res, Exception):
                if isinstance(res, FileNotFoundError) and recursive:
                    # Ignore FileNotFoundError for implicit directories returned by _expand_path.
                    if any(p.startswith(p1.rstrip("/") + "/") for p in paths1):
                        continue
                raise res

    mv = asyn.sync_wrapper(_mv)

    async def _cp_file(self, path1, path2, acl=None, **kwargs):
        """Duplicate remote file"""
        b1, k1, g1 = self.split_path(path1)
        b2, k2, g2 = self.split_path(path2)
        if g2:
            raise ValueError("Cannot write to specific object generation")
        out = await self._call(
            "POST",
            "b/{}/o/{}/rewriteTo/b/{}/o/{}",
            b1,
            k1,
            b2,
            k2,
            headers={"Content-Type": "application/json"},
            destinationPredefinedAcl=acl,
            json_out=True,
            sourceGeneration=g1,
        )
        while out["done"] is not True:
            out = await self._call(
                "POST",
                "b/{}/o/{}/rewriteTo/b/{}/o/{}",
                b1,
                k1,
                b2,
                k2,
                headers={"Content-Type": "application/json"},
                rewriteToken=out["rewriteToken"],
                destinationPredefinedAcl=acl,
                json_out=True,
                sourceGeneration=g1,
            )
        await self._write_file_cache_update(path2)

    async def _mv_file(self, path1, path2, **kwargs):
        src_bucket, src_key, generation1 = self.split_path(path1)
        dest_bucket, dest_key, generation2 = self.split_path(path2)

        if generation2:
            raise ValueError("Cannot move to specific object generation")

        if src_bucket == dest_bucket and src_key and dest_key:
            try:
                out = await self._call(
                    "POST",
                    "b/{}/o/{}/moveTo/o/{}",
                    src_bucket,
                    src_key,
                    dest_key,
                    sourceGeneration=generation1,
                    headers={
                        "Content-Type": "application/json",
                        "X-Goog-GCS-Idempotency-Token": str(uuid.uuid4()),
                    },
                    json_out=True,
                )
                await self._mv_file_cache_update(path1, path2, out)
                return
            except FileNotFoundError:
                # Raise immediately because fallback will also fail when file is not found.
                raise
            except Exception as e:
                # TODO: Fallback is added to make sure there is smooth transition, it can be removed
                # once we have metrics proving that moveTo API is working properly for all bucket types.
                logger.warning(
                    f"Failed to move file using moveTo API: {e}. Falling back to copy/delete."
                )

        await super()._mv_file(path1, path2, **kwargs)

    mv_file = asyn.sync_wrapper(_mv_file)

    async def _rm_file(self, path, **kwargs):
        bucket, key, generation = self.split_path(path)
        if key:
            await self._call("DELETE", "b/{}/o/{}", bucket, key, generation=generation)
            await self._rm_file_cache_update(path)
        else:
            await self._rmdir(path)

    async def _rm_files(self, paths):
        import random

        template = (
            "\n--===============7330845974216740156==\n"
            "Content-Type: application/http\n"
            "Content-Transfer-Encoding: binary\n"
            "Content-ID: <b29c5de2-0db4-490b-b421-6a51b598bd11+{i}>"
            "\n\nDELETE /storage/v1/b/{bucket}/o/{key}{query} HTTP/1.1\n"
            "Content-Type: application/json\n"
            "accept: application/json\ncontent-length: 0\n"
        )
        out = []
        # Splitting requests into batches
        # See https://cloud.google.com/storage/docs/batch
        for retry in range(1, 6):
            remaining = []
            chunk = paths
            parts = []
            for i, p in enumerate(chunk):
                bucket, key, generation = self.split_path(p)
                query_params = self._get_params(
                    {"generation": generation} if generation else {}
                )
                query = (
                    ("?" + "&".join(f"{k}={v}" for k, v in query_params.items()))
                    if query_params
                    else ""
                )
                parts.append(
                    template.format(
                        i=i + 1,
                        bucket=quote(bucket),
                        key=quote(key),
                        query=query,
                    )
                )
            body = "".join(parts)
            headers, content = await self._call(
                "POST",
                self.batch_url_base,
                headers={
                    "Content-Type": 'multipart/mixed; boundary="=========='
                    '=====7330845974216740156=="'
                },
                data=body + "\n--===============7330845974216740156==--",
            )

            boundary = headers["Content-Type"].split("=", 1)[1]
            txt = content.decode()
            responses = txt.split(boundary)[1:-1]
            deleted = []
            confirmed_absent = []
            for path, response in zip(paths, responses):
                m = re.search("HTTP/[0-9.]+ ([0-9]+)", response)
                code = int(m.groups()[0]) if m else None
                if code in [200, 204]:
                    out.append(path)
                    deleted.append(path)
                elif code in errs and retry < 5:
                    remaining.append(path)
                else:
                    if code == 404:
                        confirmed_absent.append(path)
                    msg = re.search("{(.*)}", response.replace("\n", ""))
                    if msg:
                        msg2 = re.search("({.*})", msg.groups()[0])
                    else:
                        msg2 = None
                    if msg and msg2:
                        out.append(OSError(msg2.groups()[0]))
                    else:
                        out.append(OSError(f"{path}: {code}"))
            # Only update the cache for objects we actually deleted, or that GCS
            # confirms are already absent. Updating it for other failed paths
            # would evict still-present objects from the cache (the HNS targeted
            # update mutates the parent listing in place and does not self-correct
            # on the next listing).
            cache_updates = deleted + confirmed_absent
            if cache_updates:
                await self._rm_files_cache_update(cache_updates)
            if remaining:
                paths = remaining
                await asyncio.sleep(min(random.random() + 2 ** (retry - 1), 32))
            else:
                break
        return out

    @property
    def on_google(self):
        # match "torage" to handle both "storage" and "Storage"
        return f"torage.{_gcp_universe_domain()}" in self._location

    async def _delete_files(self, files, batchsize):
        """Helper to delete files in batches."""
        if self.on_google:
            # emulators do not support batch
            return sum(
                await asyn._run_coros_in_chunks(
                    [
                        self._rm_files(files[i : i + batchsize])
                        for i in range(0, len(files), batchsize)
                    ],
                    return_exceptions=True,
                ),
                [],
            )
        else:
            return await asyn._run_coros_in_chunks(
                [self._rm_file(f) for f in files], return_exceptions=True, batch_size=5
            )

    async def _rm(self, path, recursive=False, maxdepth=None, batchsize=100):
        # 100 is the maximum number of operations allowed in a single GCS batch
        # request (https://cloud.google.com/storage/docs/batch); using the full
        # limit minimizes the number of round-trips when deleting many objects.
        paths = await self._expand_path(path, recursive=recursive, maxdepth=maxdepth)
        files = [p for p in paths if self.split_path(p)[1]]
        dirs = [p for p in paths if not self.split_path(p)[1]]
        exs = await self._delete_files(files, batchsize)

        # buckets
        exs.extend(
            await asyncio.gather(
                *[self._rmdir(d) for d in dirs], return_exceptions=True
            )
        )
        errors = [
            ex
            for ex in exs
            if isinstance(ex, Exception)
            and "No such object" not in str(ex)
            and not isinstance(ex, FileNotFoundError)
        ]
        if errors:
            raise errors[0]
        exs = [
            ex
            for ex in exs
            if "No such object" not in str(ex) and not isinstance(ex, FileNotFoundError)
        ]
        if not exs:
            # nothing got deleted
            raise FileNotFoundError(path)
        return exs

    rm = asyn.sync_wrapper(_rm)

    async def _pipe_file(
        self,
        path,
        data,
        metadata=None,
        consistency=None,
        content_type="application/octet-stream",
        fixed_key_metadata=None,
        chunksize=50 * 2**20,
        mode="overwrite",
    ):
        # enforce blocksize should be a multiple of 2**18
        consistency = consistency or self.consistency
        bucket, key, generation = self.split_path(path)
        size = len(data)
        out = None
        if size < chunksize:
            location = await simple_upload(
                self,
                bucket,
                key,
                data,
                metadata,
                consistency,
                content_type,
                fixed_key_metadata=fixed_key_metadata,
                mode=mode,
            )
        else:
            location = await initiate_upload(
                self,
                bucket,
                key,
                content_type,
                metadata,
                fixed_key_metadata=fixed_key_metadata,
                mode=mode,
            )
            try:
                data_view = memoryview(data)
                for offset in range(0, size, chunksize):
                    bit = data_view[offset : offset + chunksize]
                    out = await upload_chunk(
                        self, location, bit, offset, size, content_type
                    )
            except Exception:
                await self._call(
                    "DELETE",
                    location.replace("&ifGenerationMatch=0", ""),
                )
                raise

            checker = get_consistency_checker(consistency)
            checker.update(data)
            checker.validate_json_response(out)

        await self._write_file_cache_update(path)
        return location

    async def _put_file(
        self,
        lpath,
        rpath,
        metadata=None,
        consistency=None,
        content_type=None,
        chunksize=50 * 2**20,
        callback=None,
        fixed_key_metadata=None,
        mode="overwrite",
        **kwargs,
    ):
        # enforce blocksize should be a multiple of 2**18
        if os.path.isdir(lpath):
            return
        if content_type is None:
            content_type, _ = mimetypes.guess_type(lpath)
            if content_type is None:
                content_type = "application/octet-stream"
        callback = callback or NoOpCallback()
        consistency = consistency or self.consistency
        checker = get_consistency_checker(consistency)
        bucket, key, generation = self.split_path(rpath)
        if generation:
            raise ValueError("Cannot write to specific object generation")
        with open(lpath, "rb") as f0:
            size = f0.seek(0, 2)
            f0.seek(0)
            callback.set_size(size)

            if size < 5 * 2**20:
                await simple_upload(
                    self,
                    bucket,
                    key,
                    f0.read(),
                    consistency=consistency,
                    metadatain=metadata,
                    content_type=content_type,
                    fixed_key_metadata=fixed_key_metadata,
                    mode=mode,
                )
                callback.absolute_update(size)

            else:
                location = await initiate_upload(
                    self,
                    bucket,
                    key,
                    content_type,
                    metadata=metadata,
                    fixed_key_metadata=fixed_key_metadata,
                    mode=mode,
                )
                offset = 0
                try:
                    while True:
                        bit = f0.read(chunksize)
                        if not bit:
                            break
                        out = await upload_chunk(
                            self, location, bit, offset, size, content_type
                        )
                        offset += len(bit)
                        callback.absolute_update(offset)
                        checker.update(bit)
                except Exception:
                    await self._call(
                        "DELETE",
                        location.replace("&ifGenerationMatch=0", ""),
                    )
                    raise

                checker.validate_json_response(out)

            await self._write_file_cache_update(rpath)

    async def _isdir(self, path):

        try:
            return (await self._info(path))["type"] == "directory"
        except OSError:
            return False

    async def _find(
        self,
        path,
        withdirs=False,
        detail=False,
        prefix="",
        versions=False,
        maxdepth=None,
        update_cache=True,
        **kwargs,
    ):
        path = self._strip_protocol(path)

        if maxdepth is not None and maxdepth < 1:
            raise ValueError("maxdepth must be at least 1")

        # Fetch objects as if the path is a directory
        objects, _ = await self._do_list_objects(
            path, delimiter="", prefix=prefix, versions=versions
        )

        if not objects:
            # Fetch objects as if the path is a file
            bucket, key, _ = self.split_path(path)
            if prefix:
                _path = "" if not key else key.rstrip("/") + "/"
                _prefix = f"{_path}{prefix}"
            else:
                _prefix = key
            objects, _ = await self._do_list_objects(
                bucket, delimiter="", prefix=_prefix, versions=versions
            )
        else:
            _prefix = prefix

        path2 = path.rstrip("/") + "/"

        if not prefix:
            objects = [
                o for o in objects if o["name"].startswith(path2) or o["name"] == path
            ]

        dirs = self._get_dirs_and_update_cache(
            path, objects, prefix=prefix, update_cache=update_cache
        )

        if withdirs:
            objects = sorted(objects + list(dirs.values()), key=lambda x: x["name"])

        if maxdepth:
            # Filter returned objects based on requested maxdepth
            depth = path.rstrip("/").count("/") + maxdepth
            objects = list(filter(lambda o: o["name"].count("/") <= depth, objects))

        if detail:
            if versions:
                return {f"{o['name']}#{o['generation']}": o for o in objects}
            return {o["name"]: o for o in objects}

        if versions:
            return [f"{o['name']}#{o['generation']}" for o in objects]
        return [o["name"] for o in objects]

    def _get_dirs_and_update_cache(self, path, objects, prefix="", update_cache=True):
        """
        Populates the directory cache from a list of object details.

        This method reconstructs the directory hierarchy from a flat list
        of objects and update the cache, which improves the performance of
        subsequent `ls` calls.

        Parameters
        ----------
        path: str
            The root path of the find operation.
        objects: list[dict]
            A list of objects from which directories are extracted and cache is updated.
        prefix: str
            If a prefix is provided, the directory cache will *not* be updated,
            as the object list is considered partial.
        update_cache: bool
            Cache won't be updated if update_cache is False.

        Returns
        -------
        dict: A dictionary of all pseudo-directory entries created.
        """
        dirs = {}
        cache_entries = {}

        full_prefix = ""
        if prefix:
            full_prefix = path.rstrip("/") + "/" + prefix

        for obj in objects:
            # For native HNS empty folders, which are returned as directory types
            # but are not placeholders, we need to ensure they have an entry in the cache.
            if not prefix and update_cache and obj.get("type") == "directory":
                cache_entries.setdefault(obj["name"], {})

            parent = self._parent(obj["name"])
            previous = obj

            while parent:
                dir_key = self.split_path(parent)[1]
                if not dir_key or len(parent) < len(path.rstrip("/")):
                    break

                if prefix and not parent.startswith(full_prefix):
                    # If this parent doesn't match the prefix, neither will its parents.
                    break

                dirs[parent] = {
                    "Key": dir_key,
                    "Size": 0,
                    "name": parent,
                    "StorageClass": "DIRECTORY",
                    "type": "directory",
                    "size": 0,
                }

                if not prefix and update_cache:
                    listing = cache_entries.setdefault(parent, {})
                    name = previous["name"]
                    if name not in listing:
                        listing[name] = previous

                previous = dirs[parent]
                parent = self._parent(parent)
        if not prefix and update_cache:
            cache_entries_list = {k: list(v.values()) for k, v in cache_entries.items()}
            self.dircache.update(cache_entries_list)
        return dirs

    @retry_request(retries=retries)
    async def _get_file_request(
        self, rpath, lpath, *args, headers=None, callback=None, **kwargs
    ):
        rpath = self.url(rpath)
        consistency = kwargs.pop("consistency", self.consistency)
        await self._set_session()
        async with self.session.get(
            url=rpath,
            params=self._get_params(kwargs),
            headers=self._get_headers(headers),
            timeout=self.requests_timeout,
        ) as r:
            validate_response(r.status, None, rpath)
            try:
                size = int(r.headers["content-length"])
            except (KeyError, ValueError):
                size = None
            callback.set_size(size)

            checker = get_consistency_checker(consistency)
            lparent = os.path.dirname(lpath) or os.curdir
            os.makedirs(lparent, exist_ok=True)
            with open(lpath, "wb") as f2:
                while True:
                    data = await r.content.read(4096 * 32)
                    if not data:
                        break
                    f2.write(data)
                    checker.update(data)
                    callback.relative_update(len(data))

            validate_response(r.status, data, rpath)  # validate http request
            checker.validate_http_response(r)  # validate file consistency
            return r.status, r.headers, r.request_info, data

    def _init_local_file(self, lpath, total_size):
        """Creates the target directory and pre-allocates the file size."""
        os.makedirs(os.path.dirname(lpath) or os.curdir, exist_ok=True)
        if total_size == 0:
            with open(lpath, "wb"):
                pass
        else:
            with open(lpath, "wb") as f:
                f.truncate(total_size)

    @retry_request(retries=retries)
    async def _get_file_concurrent(
        self,
        rpath,
        lpath,
        concurrency,
        chunk_size,
        max_prefetch_size,
        headers=None,
        callback=None,
        fetcher_fn=None,
        **kwargs,
    ):
        """Main orchestrator for concurrent file downloads utilizing BackgroundPrefetcher."""
        details = await self._info(rpath, **kwargs)
        total_size = details.get("size", 0)

        # Concurrency typically improves performance for RAM downloads exceeding 5MB.
        # However, for disk-backed reads in the standard bucket, _cat_file uses
        # b"".join, which creates an additional data copy of chunks. Additionally,
        # prefetching is ineffective for reads under 100MB because it only activates
        # from the third read onward and scales linearly. These factors often make
        # concurrent processing slower than writing data as it arrives.
        #
        #
        # Note that the number is 5MB for zonal buckets, Thanks to our in-house, zero-copy
        # DirectMemmoveBuffer, we didn't integrated it initially with standard bucket, because
        # we first want to stabilise that in zonal bucket (lower traffic compared to standard)
        #
        # Promoting DirectMemmoveBuffer (currently used in the Zonal bucket)
        # to the standard bucket will enable in-place assembly and lower this
        # threshold. Until then, the concurrent path for standard is enabled only for disk
        # reads of 100MB or more.
        bucket, _, _ = self.split_path(rpath)
        threshold = await self._get_threshold_for_disk_reads(bucket)
        if total_size <= max(self.MIN_CHUNK_SIZE_FOR_CONCURRENCY, threshold):
            concurrency = 1

        if concurrency == 1:
            return await self._get_file_request(
                rpath, lpath, headers=headers, callback=callback, **kwargs
            )

        consistency = kwargs.pop("consistency", self.consistency)
        check_consistency = consistency not in ("none", None)

        # Prevent silent corruption by pinning the exact object generation
        generation = details.get("generation")
        if generation and "generation" not in kwargs:
            kwargs["generation"] = generation

        callback = callback or NoOpCallback()
        callback.set_size(total_size)

        # pre-allocate the file, it is required so multiple file descriptors can seek/write safely.
        self._init_local_file(lpath, total_size)
        checker = get_consistency_checker(consistency)

        if fetcher_fn is None:

            async def default_fetcher(start, size, split_factor=1):
                return await self._cat_file(
                    rpath,
                    start=start,
                    end=start + size,
                    concurrency=split_factor,
                    headers=headers,
                    **kwargs,
                )

            fetcher_fn = default_fetcher

        from .prefetcher import BackgroundPrefetcher

        prefetcher = BackgroundPrefetcher(
            fetcher=fetcher_fn,
            size=total_size,
            concurrency=concurrency,
            max_prefetch_size=max_prefetch_size,
            loop=self.loop,
        )

        fd = None

        def write_chunk(offset, chunk):
            if fd is not None:
                written = 0
                chunk_view = memoryview(chunk)
                while written < len(chunk):
                    written += os.pwrite(fd, chunk_view[written:], offset + written)
            else:
                # Thread-safe fallback for older Windows versions (Python < 3.12)
                with open(lpath, "rb+") as f:
                    f.seek(offset)
                    f.write(chunk)

        pending_writes = set()
        try:
            if hasattr(os, "pwrite"):
                fd = os.open(lpath, os.O_WRONLY | getattr(os, "O_BINARY", 0))

            async with prefetcher:
                offset = 0
                while offset < total_size:
                    read_size = min(chunk_size, total_size - offset)

                    data = await prefetcher.afetch(offset, offset + read_size)

                    if not data:
                        break

                    if check_consistency:
                        checker.update(data)

                    callback.relative_update(len(data))

                    task = asyncio.create_task(
                        asyncio.to_thread(write_chunk, offset, data)
                    )
                    pending_writes.add(task)

                    if len(pending_writes) >= concurrency:
                        done, pending_writes = await asyncio.wait(
                            pending_writes, return_when=asyncio.FIRST_COMPLETED
                        )

                        exceptions = []
                        for t in done:
                            exc = t.exception()
                            if exc:
                                exceptions.append(exc)

                        if exceptions:
                            raise exceptions[0]

                    offset += len(data)

                if offset != total_size:
                    raise aiohttp.client_exceptions.ClientError(
                        f"Expected {total_size} bytes, but only received {offset} bytes"
                    )
        finally:
            all_done = set()
            was_cancelled = False
            if pending_writes:
                while pending_writes:
                    try:
                        done_wait, pending_writes = await asyncio.wait(pending_writes)
                        all_done.update(done_wait)
                    except asyncio.CancelledError:
                        was_cancelled = True
                        pass

            if fd is not None:
                os.close(fd)

            exceptions = []
            for t in all_done:
                exc = t.exception()
                if exc:
                    exceptions.append(exc)

            if was_cancelled:
                raise asyncio.CancelledError()

            if exceptions and sys.exc_info()[1] is None:
                raise exceptions[0]

        if check_consistency:
            checker.validate_json_response(details)

    async def _get_file(self, rpath, lpath, callback=None, **kwargs):
        if os.path.isdir(lpath):
            return

        callback = callback or NoOpCallback()

        concurrency = kwargs.pop("concurrency", DEFAULT_CONCURRENCY)
        chunk_size = kwargs.pop("chunk_size", 16 * 1024 * 1024)
        max_prefetch_size = kwargs.pop(
            "max_prefetch_size", 2 * concurrency * chunk_size
        )

        try:
            # The concurrent path uses `_cat_file` to interact with gcsfs which doesn't take headers as argument.
            if concurrency > 1 and "headers" not in kwargs:
                await self._get_file_concurrent(
                    rpath,
                    lpath,
                    concurrency,
                    callback=callback,
                    chunk_size=chunk_size,
                    max_prefetch_size=max_prefetch_size,
                    **kwargs,
                )
            else:
                await self._get_file_request(rpath, lpath, callback=callback, **kwargs)
        except BaseException:
            if os.path.exists(lpath):
                os.remove(lpath)
            raise

    def _open(
        self,
        path,
        mode="rb",
        block_size=None,
        cache_options=None,
        acl=None,
        consistency=None,
        metadata=None,
        autocommit=True,
        fixed_key_metadata=None,
        generation=None,
        **kwargs,
    ):
        """
        See ``GCSFile``.

        consistency: None or str
            If None, use default for this instance
        """
        if block_size is None:
            block_size = self.default_block_size
        const = consistency or self.consistency
        return GCSFile(
            self,
            path,
            mode,
            block_size,
            cache_options=cache_options,
            consistency=const,
            metadata=metadata,
            acl=acl,
            autocommit=autocommit,
            fixed_key_metadata=fixed_key_metadata,
            generation=generation,
            **kwargs,
        )

    @classmethod
    def _split_path(cls, path, version_aware=False):
        """
        Normalise GCS path string into bucket and key.

        Parameters
        ----------
        path : string
            Input path, like `gcs://mybucket/path/to/file`.
            Path is of the form: '[gs|gcs://]bucket[/key][?querystring][#fragment]'

        GCS allows object generation (object version) to be specified in either
        the URL fragment or the `generation` query parameter. When provided,
        the fragment will take priority over the `generation` query parameter.

        Returns
        -------
            (bucket, key, generation) tuple
        """
        path = cls._strip_protocol(path).lstrip("/")
        if "/" not in path:
            return path, "", None
        bucket, keypart = path.split("/", 1)
        key = keypart
        generation = None
        if version_aware:
            parts = urlsplit(keypart)
            try:
                if parts.fragment:
                    generation = parts.fragment
                elif parts.query:
                    parsed = parse_qs(parts.query)
                    if "generation" in parsed:
                        generation = parsed["generation"][0]
                # Sanity check whether this could be a valid generation ID. If
                # it is not, assume that # or ? characters are supposed to be
                # part of the object name.
                if generation is not None:
                    int(generation)
                    key = parts.path
            except ValueError:
                generation = None
        return (
            bucket,
            key,
            generation,
        )

    def split_path(self, path):
        return self._split_path(path, version_aware=self.version_aware)

    def sign(self, path, expiration=100, **kwargs):
        """Create a signed URL representing the given path.

        Parameters
        ----------
        path : str
             The path on the filesystem
        expiration : int
            Number of seconds to enable the URL for

        Returns
        -------
        URL : str
            The signed URL
        """
        from google.cloud import storage

        client = storage.Client(
            credentials=self.credentials.credentials,
            project=self.project,
        )

        bucket, key, generation = self.split_path(path)
        bucket = client.bucket(bucket)
        blob = bucket.blob(key)

        return blob.generate_signed_url(
            expiration=timedelta(seconds=expiration),
            generation=generation,
            api_access_endpoint=self._endpoint,
            **kwargs,
        )


GoogleCredentials.load_tokens()


class GCSFile(fsspec.spec.AbstractBufferedFile):
    def __init__(
        self,
        gcsfs,
        path,
        mode="rb",
        block_size=DEFAULT_BLOCK_SIZE,
        autocommit=True,
        cache_type="readahead",
        cache_options=None,
        acl=None,
        consistency="md5",
        metadata=None,
        content_type=None,
        timeout=None,
        fixed_key_metadata=None,
        generation=None,
        kms_key_name=None,
        **kwargs,
    ):
        """
        Open a file.

        Parameters
        ----------
        gcsfs: instance of GCSFileSystem
        path: str
            location in GCS, like 'bucket/path/to/file'
        mode: str
            Normal file modes. Currently only 'wb' amd 'rb'.
        block_size: int
            Buffer size for reading or writing
        acl: str
            ACL to apply, if any, one of ``ACLs``. New files are normally
            "bucketownerfullcontrol", but a default can be configured per
            bucket.
        consistency: str, 'none', 'size', 'md5', 'crc32c'
            Check for success in writing, applied at file close.
            'size' ensures that the number of bytes reported by GCS matches
            the number we wrote; 'md5' does a full checksum. Any value other
            than 'size' or 'md5' or 'crc32c' is assumed to mean no checking.
        content_type: str
            default when unspecified is provided by mimetypes.guess_type or
            otherwise `application/octet-stream`. See the list of available
            content types at https://www.iana.org/assignments/media-types/media-types.txt
        metadata: dict
            Custom metadata, in key/value pairs, added at file creation
        fixed_key_metadata: dict
            Google metadata, in key/value pairs, supported keys:
                - cache_control
                - content_disposition
                - content_encoding
                - content_language
                - custom_time
            More info:
            https://cloud.google.com/storage/docs/metadata#mutable
        kms_key_name: str
            Resource name of the Cloud KMS key that will be used to encrypt
            the object.
            More info:
            https://cloud.google.com/storage/docs/encryption/customer-managed-keys
        timeout: int
            Timeout seconds for the asynchronous callback.
        generation: str
            Object generation.
        """
        bucket, key, path_generation = gcsfs.split_path(path)
        if not key:
            raise OSError("Attempt to open a bucket")
        self.generation = _coalesce_generation(generation, path_generation)
        self.concurrency = kwargs.get("concurrency", DEFAULT_CONCURRENCY)
        super().__init__(
            gcsfs,
            path,
            mode,
            block_size,
            autocommit=autocommit,
            cache_type=cache_type,
            cache_options=cache_options,
            **kwargs,
        )
        self.cache_type = cache_type
        self.gcsfs = gcsfs
        self.bucket = bucket
        self.key = key
        self.acl = acl
        self.consistency = consistency
        self.checker = get_consistency_checker(consistency)

        # Ideally, all of these fields should be part of `cache_options`. Because current
        # `fsspec` caches do not accept arbitrary `*args` and `**kwargs`, passing them
        # there currently causes instantiation errors. We are holding off on introducing
        # them as explicit keyword arguments to ensure existing user workloads are not
        # disrupted. This will be refactored once the upstream `fsspec` changes are merged.
        use_prefetch_reader = kwargs.get(
            "use_experimental_adaptive_prefetching", False
        ) or os.environ.get(
            "USE_EXPERIMENTAL_ADAPTIVE_PREFETCHING", "false"
        ).lower() in (
            "true",
            "1",
        )

        if "r" in mode and use_prefetch_reader:
            max_prefetch_size = kwargs.get("max_prefetch_size", MAX_PREFETCH_SIZE)
            from .prefetcher import BackgroundPrefetcher

            self._prefetch_engine = BackgroundPrefetcher(
                self._async_fetch_range,
                self.size,
                max_prefetch_size=max_prefetch_size,
                concurrency=self.concurrency,
                loop=self.gcsfs.loop,
            )
        else:
            self._prefetch_engine = None

        # _supports_append is an internal argument not meant to be used directly.
        # If True, allows opening file in append mode. This is generally not supported
        # by GCS, but may be supported by subclasses (e.g. ZonalFile). This flag should
        # be set by subclasses that support append operations. Otherwise, the mode
        # will be overwritten to "wb" mode with a warning.
        _supports_append = kwargs.pop("_supports_append", False)
        if "a" in self.mode and not _supports_append:
            warnings.warn(
                "Append mode 'a' is not supported in GCS. Using overwrite mode instead."
            )
            self.mode = self.mode.replace("a", "w")

        if "r" in self.mode:
            det = self.details
        else:
            det = {}
        self.content_type = content_type or det.get(
            "contentType",
            mimetypes.guess_type(self.path)[0] or "application/octet-stream",
        )
        self.metadata = metadata or det.get("metadata", {})
        self.fixed_key_metadata = _convert_fixed_key_metadata(det, from_google=True)
        self.fixed_key_metadata.update(fixed_key_metadata or {})
        self.kms_key_name = kms_key_name
        self.timeout = timeout
        if mode in {"wb", "xb"}:
            if self.blocksize < GCS_MIN_BLOCK_SIZE:
                warnings.warn("Setting block size to minimum value, 2**18")
                self.blocksize = GCS_MIN_BLOCK_SIZE
            self.location = None

    @property
    def details(self):
        if self._details is None:
            self._details = self.fs.info(self.path, generation=self.generation)
        return self._details

    def info(self):
        """File information about this path"""
        return self.details

    def url(self):
        """HTTP link to this file's data"""
        return self.fs.url(self.path, generation=self.generation)

    def _upload_chunk(self, final=False):
        """Write one part of a multi-block file upload

        Parameters
        ----------
        final: bool
            Complete and commit upload
        """
        while True:
            # shortfall splits blocks bigger than max allowed upload
            data_view = self.buffer.getbuffer()
            head = {}
            l = len(data_view)

            if (l < GCS_MIN_BLOCK_SIZE) and (not final or not self.autocommit):
                # either flush() was called, but we don't have enough to
                # push, or we split a big upload, and have less left than one
                # block.  If this is the final part, OK to violate those
                # terms.
                return False

            # Select the biggest possible chunk of data to be uploaded
            chunk_length = min(l, GCS_MAX_BLOCK_SIZE)
            # This chunk finalizes the upload when it is the final flush, we are
            # autocommitting, and all remaining data fits in a single chunk.
            finalizes_upload = final and self.autocommit and chunk_length == l
            if not finalizes_upload:
                # GCS requires non-final resumable-upload chunks to be
                # multiples of 256 KiB:
                # https://cloud.google.com/storage/docs/performing-resumable-uploads#multiple-chunk-upload
                chunk_length = (chunk_length // GCS_MIN_BLOCK_SIZE) * GCS_MIN_BLOCK_SIZE
            chunk = data_view[:chunk_length]
            if finalizes_upload:
                if l:
                    # last chunk
                    head["Content-Range"] = "bytes %i-%i/%i" % (
                        self.offset,
                        self.offset + chunk_length - 1,
                        self.offset + l,
                    )
                else:
                    # closing when buffer is empty
                    head["Content-Range"] = "bytes */%i" % self.offset
                    chunk = None
            else:
                head["Content-Range"] = "bytes %i-%i/*" % (
                    self.offset,
                    self.offset + chunk_length - 1,
                )
            head.update(
                {"Content-Type": self.content_type, "Content-Length": str(chunk_length)}
            )
            headers, contents = self.gcsfs.call(
                "POST", self.location, headers=head, data=chunk
            )
            if "Range" in headers:
                end = int(headers["Range"].split("-")[1])
                shortfall = (self.offset + l - 1) - end
                if shortfall > 0:
                    self.checker.update(data_view[:-shortfall])
                    self.buffer = UnclosableBytesIO(data_view[-shortfall:])
                    self.buffer.seek(shortfall)
                    self.offset += l - shortfall
                    continue
                else:
                    self.checker.update(data_view)
                if final and contents:
                    j = json.loads(contents)
                    self.generation = j.get("generation")
            else:
                assert final, "Response looks like upload is over"
                if l:
                    j = json.loads(contents)
                    self.checker.update(data_view)
                    self.checker.validate_json_response(j)
                    self.generation = j.get("generation")
            # Clear buffer and update offset when all is received
            self.buffer = UnclosableBytesIO()
            self.offset += l
            break
        return True

    def commit(self):
        """If not auto-committing, finalize file"""
        self.autocommit = True
        self._upload_chunk(final=True)

    def _initiate_upload(self):
        """Create multi-upload"""
        self.location = asyn.sync(
            self.gcsfs.loop,
            initiate_upload,
            self.gcsfs,
            self.bucket,
            self.key,
            self.content_type,
            self.metadata,
            self.fixed_key_metadata,
            mode="create" if "x" in self.mode else "overwrite",
            kms_key_name=self.kms_key_name,
            timeout=self.timeout,
        )

    def discard(self):
        """Cancel in-progress multi-upload

        Should only happen during discarding this write-mode file
        """
        if self.location is None:
            return
        self.gcsfs.call(
            "DELETE",
            self.location.replace("&ifGenerationMatch=0", ""),
        )
        self.location = None
        self.closed = True

    def _simple_upload(self):
        """One-shot upload, less than 5MB"""
        self.buffer.seek(0)
        data = self.buffer.read()
        j = asyn.sync(
            self.gcsfs.loop,
            simple_upload,
            self.gcsfs,
            self.bucket,
            self.key,
            data,
            self.metadata,
            self.consistency,
            self.content_type,
            self.fixed_key_metadata,
            mode="create" if "x" in self.mode else "overwrite",
            kms_key_name=self.kms_key_name,
            timeout=self.timeout,
        )
        self.generation = j.get("generation")

    def _fetch_range(self, start=None, end=None):
        """Get data from GCS

        start, end : None or integers
            if not both None, fetch only given range
        """
        try:
            if hasattr(self, "_prefetch_engine") and self._prefetch_engine:
                return self._prefetch_engine.fetch(start=start, end=end)
            return self.fs.cat_file(
                self.path,
                start=start,
                end=end,
                concurrency=self.concurrency,
                cache_type=self.cache_type,
            )
        except RuntimeError as e:
            if "not satisfiable" in str(e):
                return b""
            raise

    async def _async_fetch_range(self, start_offset, total_size, split_factor=1):
        """Async fetcher mapped to the Prefetcher engine for regional buckets."""
        return await self.gcsfs._cat_file_concurrent(
            self.path,
            start=start_offset,
            end=start_offset + total_size,
            concurrency=split_factor,
            cache_type=self.cache_type,
        )

    def close(self):
        super().close()
        if hasattr(self, "_prefetch_engine") and self._prefetch_engine:
            self._prefetch_engine.close()


def _convert_fixed_key_metadata(metadata, *, from_google=False):
    """
    Convert fixed_key_metadata to/from GCS format.

    Google uses camelCase for its parameters, this function transform
    exposed fixed_key_metadata (snake_case) to or from GCS(google) format

    Parameters
    ----------
    metadata: dict
        A key value pair of fixed_key_metadata, key can be either
        camel case or snake case.
    from_google: bool
        True means that the metadata come from google and thus should be converted
        to snake_case
    """
    out = {}
    if metadata is None:
        return out

    for key, attribute_name in SUPPORTED_FIXED_KEY_METADATA.items():
        src = key if not from_google else attribute_name
        dst = attribute_name if not from_google else key
        if src in metadata:
            out[dst] = metadata[src]
    return out


async def upload_chunk(fs, location, data, offset, size, content_type):
    """
    Uploads a chunk of data. This function has a conditional path to support
    experimental features for Zonal buckets to append data using gRPC.
    """
    from google.cloud.storage.asyncio.async_appendable_object_writer import (
        AsyncAppendableObjectWriter,
    )

    from .extended_gcsfs import ExtendedGcsFileSystem
    from .extended_gcsfs import upload_chunk as ext_upload_chunk

    # location is AsyncAppendableObjectWriter only when ExtendedGcsFileSystem is used
    if isinstance(fs, ExtendedGcsFileSystem) and isinstance(
        location, AsyncAppendableObjectWriter
    ):

        return await ext_upload_chunk(fs, location, data, offset, size, content_type)
    head = {}
    l = len(data)
    range = "bytes %i-%i/%i" % (offset, offset + l - 1, size)
    head["Content-Range"] = range
    head.update({"Content-Type": content_type, "Content-Length": str(l)})
    # aiohttp handles bytes and memoryview natively and zero-copy;
    # no need to wrap in UnclosableBytesIO.
    payload = (
        data
        if isinstance(data, (bytes, bytearray, memoryview))
        else UnclosableBytesIO(data)
    )
    headers, txt = await fs._call("POST", location, headers=head, data=payload)
    if "Range" in headers:
        end = int(headers["Range"].split("-")[1])
        shortfall = (offset + l - 1) - end
        if shortfall:
            return await upload_chunk(
                fs, location, data[-shortfall:], end + 1, size, content_type
            )
    return json.loads(txt) if txt else None


async def initiate_upload(
    fs,
    bucket,
    key,
    content_type="application/octet-stream",
    metadata=None,
    fixed_key_metadata=None,
    mode="overwrite",
    kms_key_name=None,
):
    """
    Initiates a resumable upload. This function has a conditional path to support
    experimental features for Zonal buckets to append data using gRPC, returning an
    "AsyncAppendableObjectWriter" instance as location.
    """
    from .extended_gcsfs import ExtendedGcsFileSystem
    from .extended_gcsfs import initiate_upload as ext_initiate_upload

    # Explicit type checking is used to ensure only the ExtendedGcsFileSystem
    # enters this path, ruling out false positives from mocks or coincidentally matching attributes.
    if isinstance(fs, ExtendedGcsFileSystem) and await fs._is_zonal_bucket(bucket):

        return await ext_initiate_upload(
            fs,
            bucket,
            key,
            content_type,
            metadata,
            fixed_key_metadata,
            mode,
            kms_key_name,
        )

    j = {"name": key}
    if metadata:
        j["metadata"] = metadata
    kw = {"ifGenerationMatch": "0"} if mode == "create" else {}
    if kms_key_name:
        kw["kmsKeyName"] = kms_key_name
    j.update(_convert_fixed_key_metadata(fixed_key_metadata))
    headers, _ = await fs._call(
        method="POST",
        path=f"{fs._location}/upload/storage/v1/b/{quote(bucket)}/o?name={quote(key)}",
        uploadType="resumable",
        json=j,
        headers={"X-Upload-Content-Type": content_type},
        **kw,
    )
    loc = headers["Location"]
    out = loc[0] if isinstance(loc, list) else loc  # <- for CVR responses
    if len(str(loc)) < 20:
        logger.error("Location failed: %s" % headers)
    return out


async def simple_upload(
    fs,
    bucket,
    key,
    datain,
    metadatain=None,
    consistency=None,
    content_type="application/octet-stream",
    fixed_key_metadata=None,
    mode="overwrite",
    kms_key_name=None,
):
    """
    Performs a simple, single-request upload. This function has a conditional path to support
    experimental features for Zonal buckets to upload data using gRPC.
    """
    from .extended_gcsfs import ExtendedGcsFileSystem
    from .extended_gcsfs import simple_upload as ext_simple_upload

    # Explicit type checking is used to ensure only the ExtendedGcsFileSystem
    # enters this path, ruling out false positives from mocks or coincidentally matching attributes.
    if isinstance(fs, ExtendedGcsFileSystem) and await fs._is_zonal_bucket(bucket):

        return await ext_simple_upload(
            fs,
            bucket,
            key,
            datain,
            metadatain,
            consistency,
            content_type,
            fixed_key_metadata,
            mode,
            kms_key_name,
        )

    checker = get_consistency_checker(consistency)
    path = f"{fs._location}/upload/storage/v1/b/{quote(bucket)}/o"
    metadata = {"name": key}
    if metadatain is not None:
        metadata["metadata"] = metadatain
    kw = {"ifGenerationMatch": "0"} if mode == "create" else {}
    if kms_key_name:
        kw["kmsKeyName"] = kms_key_name
    metadata.update(_convert_fixed_key_metadata(fixed_key_metadata))
    metadata = json.dumps(metadata)
    template = (
        "--==0=="
        "\nContent-Type: application/json; charset=UTF-8"
        "\n\n" + metadata + "\n--==0==" + f"\nContent-Type: {content_type}" + "\n\n"
    )

    data = template.encode() + datain + b"\n--==0==--"
    # aiohttp handles bytes and memoryview natively and zero-copy;
    # no need to wrap in UnclosableBytesIO.
    payload = (
        data
        if isinstance(data, (bytes, bytearray, memoryview))
        else UnclosableBytesIO(data)
    )
    j = await fs._call(
        "POST",
        path,
        uploadType="multipart",
        headers={"Content-Type": 'multipart/related; boundary="==0=="'},
        data=payload,
        json_out=True,
        **kw,
    )
    checker.update(datain)
    checker.validate_json_response(j)
    return j

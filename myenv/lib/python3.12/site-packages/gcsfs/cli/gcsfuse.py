import logging

import click
from fuse import FUSE

from gcsfs.gcsfuse import GCSFS


@click.command()
@click.argument("bucket", type=str, required=True)
@click.argument("mount_point", type=str, required=True)
@click.option(
    "--token",
    type=str,
    required=False,
    default=None,
    help="Token to use for authentication",
)
@click.option(
    "--project-id", type=str, required=False, default="", help="Billing Project ID"
)
@click.option(
    "--foreground/--background",
    default=True,
    help="Run in the foreground or as a background process",
)
@click.option(
    "--threads/--no-threads", default=True, help="Whether to run with threads"
)
@click.option(
    "--cache_files", type=int, default=10, help="Number of open files to cache"
)
@click.option(
    "-v",
    "--verbose",
    count=True,
    help="Set logging level. '-v' for 'gcsfuse' logging."
    "'-v -v' for complete debug logging.",
)
def main(
    bucket, mount_point, token, project_id, foreground, threads, cache_files, verbose
):
    """Mount a Google Cloud Storage (GCS) bucket to a local directory"""

    if verbose == 1:
        logging.basicConfig(level=logging.INFO)
        logging.getLogger("gcsfs.gcsfuse").setLevel(logging.DEBUG)
    if verbose > 1:
        logging.basicConfig(level=logging.DEBUG)

    fmt = "%(asctime)s %(name)-12s %(levelname)-8s %(message)s"
    if verbose == 1:
        logging.basicConfig(level=logging.INFO, format=fmt)
        logging.getLogger("gcsfs.gcsfuse").setLevel(logging.DEBUG)
    if verbose > 1:
        logging.basicConfig(level=logging.DEBUG, format=fmt)

    print(f"Mounting bucket {bucket} to directory {mount_point}")
    print("foreground:", foreground, ", nothreads:", not threads)
    FUSE(
        GCSFS(bucket, token=token, project=project_id, nfiles=cache_files),
        mount_point,
        nothreads=not threads,
        foreground=foreground,
    )


if __name__ == "__main__":
    main()

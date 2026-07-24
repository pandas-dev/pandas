"""The cli for auth."""

import argparse
import sys
import warnings
from getpass import getpass

from jupyter_core.paths import jupyter_config_dir
from traitlets.log import get_logger

from jupyter_server.auth import passwd  # type:ignore[attr-defined]
from jupyter_server.config_manager import BaseJSONConfigManager


def set_password(args):
    """Set a password."""
    password = args.password

    while not password:
        password1 = getpass("" if args.quiet else "Provide password: ")
        password_repeat = getpass("" if args.quiet else "Repeat password:  ")
        if password1 != password_repeat:
            warnings.warn("Passwords do not match, try again", stacklevel=2)
        elif len(password1) < 4:
            warnings.warn("Please provide at least 4 characters", stacklevel=2)
        else:
            password = password1

    password_hash = passwd(password)
    cfg = BaseJSONConfigManager(config_dir=jupyter_config_dir())
    cfg.update(
        "jupyter_server_config",
        {
            "ServerApp": {
                "password": password_hash,
            }
        },
    )
    if not args.quiet:
        log = get_logger()
        log.info("password stored in config dir: %s" % jupyter_config_dir())


def main(argv):
    """The main cli handler."""
    parser = argparse.ArgumentParser(argv[0])
    subparsers = parser.add_subparsers()
    parser_password = subparsers.add_parser(
        "password", help="sets a password for your jupyter server"
    )
    parser_password.add_argument(
        "password",
        help="password to set, if not given, a password will be queried for (NOTE: this may not be safe)",
        nargs="?",
    )
    parser_password.add_argument("--quiet", help="suppress messages", action="store_true")
    parser_password.set_defaults(function=set_password)
    args = parser.parse_args(argv[1:])
    args.function(args)


if __name__ == "__main__":
    main(sys.argv)

import argparse
import logging
import os
import signal
import sys
from http.server import ThreadingHTTPServer
from typing import Any

from moto.moto_proxy import logger
from moto.moto_proxy.proxy3 import CertificateCreator, ProxyRequestHandler, with_color


def signal_handler(signum: Any, frame: Any) -> None:  # pylint: disable=unused-argument
    sys.exit(0)


def get_help_msg() -> str:
    msg = """
    ###################################################################################
    $$___$$_ __$$$___ $$$$$$_ __$$$___\t__$$$$$$__ $$$$$$__ __$$$___ $$___$$_ $$____$$_
    $$$_$$$_ _$$_$$__ __$$___ _$$_$$__\t__$$___$$_ $$___$$_ _$$_$$__ $$$_$$$_ _$$__$$__
    $$$$$$$_ $$___$$_ __$$___ $$___$$_\t__$$___$$_ $$___$$_ $$___$$_ _$$$$$__ __$$$$___
    $$_$_$$_ $$___$$_ __$$___ $$___$$_\t__$$$$$$__ $$$$$$__ $$___$$_ _$$$$$__ ___$$____
    $$___$$_ _$$_$$__ __$$___ _$$_$$__\t__$$______ $$___$$_ _$$_$$__ $$$_$$$_ ___$$____
    $$___$$_ __$$$___ __$$___ __$$$___\t__$$______ $$___$$_ __$$$___ $$___$$_ ___$$____
    ###################################################################################"""
    msg += "\n"
    msg += "Using the CLI:"
    msg += "\n"
    msg += with_color(37, text="\texport HTTPS_PROXY=http://localhost:5005")
    msg += "\n"
    msg += with_color(37, text="\taws cloudformation list-stacks --no-verify-ssl\n")
    msg += "\n"
    msg += "Using pytest:"
    msg += "\n"
    msg += with_color(37, text=f"\texport AWS_CA_BUNDLE={CertificateCreator.cacert}")
    msg += "\n"
    msg += with_color(
        37,
        text="\tHTTPS_PROXY=http://localhost:5005 MOTO_PROXY_PORT=5005 pytest tests_dir\n",
    )
    return msg


def main(argv: Any = None) -> None:
    argv = argv or sys.argv[1:]
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter, description=get_help_msg()
    )

    parser.add_argument(
        "-H", "--host", type=str, help="Which host to bind", default="127.0.0.1"
    )
    parser.add_argument(
        "-p",
        "--port",
        type=int,
        help="Port number to use for connection",
        default=int(os.environ.get("MOTO_PROXY_PORT", 5005)),
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Add verbose logging",
    )

    args = parser.parse_args(argv)

    if args.verbose:
        logger.setLevel(logging.DEBUG)

    ProxyRequestHandler.validate()

    if "MOTO_PORT" not in os.environ:
        os.environ["MOTO_PORT"] = f"{args.port}"
    os.environ["TEST_PROXY_MODE"] = "true"

    try:
        signal.signal(signal.SIGINT, signal_handler)
        signal.signal(signal.SIGTERM, signal_handler)
    except Exception:
        pass  # ignore "ValueError: signal only works in main thread"

    server_address = (args.host, args.port)

    httpd = ThreadingHTTPServer(server_address, ProxyRequestHandler)

    sa = httpd.socket.getsockname()

    print("Call `moto_proxy -h` for example invocations")  # noqa: T201
    print(f"Serving HTTP Proxy on {sa[0]}:{sa[1]} ...")  # noqa: T201
    httpd.serve_forever()


if __name__ == "__main__":
    main()

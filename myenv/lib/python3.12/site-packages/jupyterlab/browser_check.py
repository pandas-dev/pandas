# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.

"""
This module is meant to run JupyterLab in a headless browser, making sure
the application launches and starts up without errors.
"""

import asyncio
import inspect
import logging
import os
import shutil
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from os import path as osp

from jupyter_server.serverapp import aliases, flags
from jupyter_server.utils import pathname2url, urljoin
from tornado.ioloop import IOLoop
from tornado.iostream import StreamClosedError
from tornado.websocket import WebSocketClosedError
from traitlets import Bool, Unicode

from .labapp import LabApp, get_app_dir
from .tests.test_app import TestEnv

here = osp.abspath(osp.dirname(__file__))
test_flags = dict(flags)
test_flags["core-mode"] = ({"BrowserApp": {"core_mode": True}}, "Start the app in core mode.")
test_flags["dev-mode"] = ({"BrowserApp": {"dev_mode": True}}, "Start the app in dev mode.")
test_flags["watch"] = ({"BrowserApp": {"watch": True}}, "Start the app in watch mode.")

test_aliases = dict(aliases)
test_aliases["app-dir"] = "BrowserApp.app_dir"


class LogErrorHandler(logging.Handler):
    """A handler that exits with 1 on a logged error."""

    def __init__(self):
        super().__init__(level=logging.ERROR)
        self.errored = False

    def filter(self, record):
        # Handle known StreamClosedError from Tornado
        # These occur when we forcibly close Websockets or
        # browser connections during the test.
        # https://github.com/tornadoweb/tornado/issues/2834
        if (
            hasattr(record, "exc_info")
            and record.exc_info is not None
            and isinstance(record.exc_info[1], (StreamClosedError, WebSocketClosedError))
        ):
            return
        return super().filter(record)

    def emit(self, record):
        print(record.msg, file=sys.stderr)
        self.errored = True


def run_test(app, func):
    """Synchronous entry point to run a test function.
    func is a function that accepts an app url as a parameter and returns a result.
    func can be synchronous or asynchronous.  If it is synchronous, it will be run
    in a thread, so asynchronous is preferred.
    """
    IOLoop.current().spawn_callback(run_test_async, app, func)


async def run_test_async(app, func):
    """Run a test against the application.
    func is a function that accepts an app url as a parameter and returns a result.
    func can be synchronous or asynchronous.  If it is synchronous, it will be run
    in a thread, so asynchronous is preferred.
    """
    handler = LogErrorHandler()
    app.log.addHandler(handler)

    env_patch = TestEnv()
    env_patch.start()

    app.log.info("Running async test")

    # The entry URL for browser tests is different in notebook >= 6.0,
    # since that uses a local HTML file to point the user at the app.
    if hasattr(app, "browser_open_file"):
        url = urljoin("file:", pathname2url(app.browser_open_file))
    else:
        url = app.display_url

    # Allow a synchronous function to be passed in.
    if inspect.iscoroutinefunction(func):
        test = func(url)
    else:
        app.log.info("Using thread pool executor to run test")
        loop = asyncio.get_event_loop()
        executor = ThreadPoolExecutor()
        task = loop.run_in_executor(executor, func, url)
        test = asyncio.wait([task])

    try:
        await test
    except Exception as e:
        app.log.critical("Caught exception during the test:")
        app.log.error(str(e))

    app.log.info("Test Complete")

    result = 0
    if handler.errored:
        result = 1
        app.log.critical("Exiting with 1 due to errors")
    else:
        app.log.info("Exiting normally")

    app.log.info("Stopping server...")
    try:
        app.http_server.stop()
        app.io_loop.stop()
        env_patch.stop()
    except Exception as e:
        app.log.error(str(e))
        result = 1
    finally:
        time.sleep(2)
        os._exit(result)


async def run_async_process(cmd, **kwargs):
    """Run an asynchronous command"""
    proc = await asyncio.create_subprocess_exec(*cmd, **kwargs)
    stdout, stderr = await proc.communicate()
    if proc.returncode != 0:
        raise RuntimeError(str(cmd) + " exited with " + str(proc.returncode))
    return stdout, stderr


async def run_browser(url):
    """Run the browser test and return an exit code."""
    target = osp.join(get_app_dir(), "browser_test")
    if not osp.exists(osp.join(target, "node_modules")):
        if not osp.exists(target):
            os.makedirs(osp.join(target))
        await run_async_process(["npm", "init", "-y"], cwd=target)
        await run_async_process(["npm", "install", "playwright@^1.9.2"], cwd=target)
    await run_async_process(["npx", "playwright", "install"], cwd=target)
    shutil.copy(osp.join(here, "browser-test.js"), osp.join(target, "browser-test.js"))
    await run_async_process(["node", "browser-test.js", url], cwd=target)


def run_browser_sync(url):
    """Run the browser test and return an exit code."""
    target = osp.join(get_app_dir(), "browser_test")
    if not osp.exists(osp.join(target, "node_modules")):
        os.makedirs(target)
        subprocess.call(["npm", "init", "-y"], cwd=target)  # noqa S603 S607
        subprocess.call(["npm", "install", "playwright@^1.9.2"], cwd=target)  # noqa S603 S607
    subprocess.call(["npx", "playwright", "install"], cwd=target)  # noqa S603 S607
    shutil.copy(osp.join(here, "browser-test.js"), osp.join(target, "browser-test.js"))
    return subprocess.check_call(["node", "browser-test.js", url], cwd=target)  # noqa S603 S607


class BrowserApp(LabApp):
    """An app the launches JupyterLab and waits for it to start up, checking for
    JS console errors, JS errors, and Python logged errors.
    """

    name = __name__
    open_browser = False

    serverapp_config = {"base_url": "/foo/"}
    default_url = Unicode("/lab?reset", config=True, help="The default URL to redirect to from `/`")
    ip = "127.0.0.1"
    flags = test_flags
    aliases = test_aliases
    test_browser = Bool(True)

    def initialize_settings(self):
        self.settings.setdefault("page_config_data", {})
        self.settings["page_config_data"]["browserTest"] = True
        self.settings["page_config_data"]["buildAvailable"] = False
        self.settings["page_config_data"]["exposeAppInBrowser"] = True
        super().initialize_settings()

    def initialize_handlers(self):
        func = run_browser if self.test_browser else lambda url: 0
        if os.name == "nt" and func == run_browser:
            func = run_browser_sync
        run_test(self.serverapp, func)
        super().initialize_handlers()


def _jupyter_server_extension_points():
    return [{"module": __name__, "app": BrowserApp}]


def _jupyter_server_extension_paths():
    return [{"module": "jupyterlab.browser_check"}]


if __name__ == "__main__":
    skip_options = ["--no-browser-test", "--no-chrome-test"]
    for option in skip_options:
        if option in sys.argv:
            BrowserApp.test_browser = False
            sys.argv.remove(option)

    BrowserApp.launch_instance()

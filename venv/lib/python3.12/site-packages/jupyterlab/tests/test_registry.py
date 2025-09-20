"""Test yarn registry replacement"""

# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.

import logging
import subprocess
from os.path import join as pjoin
from unittest.mock import patch

from jupyterlab import commands

from .test_jupyterlab import AppHandlerTest


class TestAppHandlerRegistry(AppHandlerTest):
    def test_node_not_available(self):
        # patch should be applied on `jupyterlab.commands` and not on `jupyterlab_server.process`
        # See https://docs.python.org/3/library/unittest.mock.html#where-to-patch
        with patch("jupyterlab.commands.which") as which:
            which.side_effect = ValueError("Command not found")

            logger = logging.getLogger("jupyterlab")
            config = commands._yarn_config(logger)

            which.assert_called_once_with("node")
            self.assertDictEqual(config, {"yarn config": {}, "npm config": {}})

    def test_yarn_config(self):
        with patch("subprocess.check_output") as check_output:
            yarn_registry = "https://private.yarn/manager"
            check_output.return_value = b"\n".join(
                [
                    b'{"type":"info","data":"yarn config"}',
                    b'{"type":"inspect","data":{"registry":"'
                    + bytes(yarn_registry, "utf-8")
                    + b'"}}',
                    b'{"type":"info","data":"npm config"}',
                    b'{"type":"inspect","data":{"registry":"'
                    + bytes(yarn_registry, "utf-8")
                    + b'"}}',
                ]
            )
            logger = logging.getLogger("jupyterlab")
            config = commands._yarn_config(logger)

            self.assertDictEqual(
                config,
                {
                    "yarn config": {"registry": yarn_registry},
                    "npm config": {"registry": yarn_registry},
                },
            )

    def test_yarn_config_failure(self):
        with patch("subprocess.check_output") as check_output:
            check_output.side_effect = subprocess.CalledProcessError(
                1, ["yarn", "config", "list"], b"", stderr=b"yarn config failed."
            )

            logger = logging.getLogger("jupyterlab")
            config = commands._yarn_config(logger)

            self.assertDictEqual(config, {"yarn config": {}, "npm config": {}})

    def test_get_registry(self):
        with patch("subprocess.check_output") as check_output:
            yarn_registry = "https://private.yarn/manager"
            check_output.return_value = b"\n".join(
                [
                    b'{"type":"info","data":"yarn config"}',
                    b'{"type":"inspect","data":{"registry":"'
                    + bytes(yarn_registry, "utf-8")
                    + b'"}}',
                    b'{"type":"info","data":"npm config"}',
                    b'{"type":"inspect","data":{"registry":"'
                    + bytes(yarn_registry, "utf-8")
                    + b'"}}',
                ]
            )

            handler = commands.AppOptions()

            self.assertEqual(handler.registry, yarn_registry)

    def test_populate_staging(self):
        with patch("subprocess.check_output") as check_output:
            yarn_registry = "https://private.yarn/manager"
            check_output.return_value = b"\n".join(
                [
                    b'{"type":"info","data":"yarn config"}',
                    b'{"type":"inspect","data":{"registry":"'
                    + bytes(yarn_registry, "utf-8")
                    + b'"}}',
                    b'{"type":"info","data":"npm config"}',
                    b'{"type":"inspect","data":{"registry":"'
                    + bytes(yarn_registry, "utf-8")
                    + b'"}}',
                ]
            )

            staging = pjoin(self.app_dir, "staging")
            handler = commands._AppHandler(commands.AppOptions())
            handler._populate_staging()

            lock_path = pjoin(staging, "yarn.lock")
            with open(lock_path) as f:
                lock = f.read()

            # yarn >=2.x does not record the registry in the lockfile
            self.assertNotIn(commands.YARN_DEFAULT_REGISTRY, lock)
            self.assertNotIn(yarn_registry, lock)

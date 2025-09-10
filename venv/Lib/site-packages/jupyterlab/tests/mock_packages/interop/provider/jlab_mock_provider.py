# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.

import json
import os.path as osp

HERE = osp.abspath(osp.dirname(__file__))

with open(osp.join(HERE, "static", "package.json")) as fid:
    data = json.load(fid)


def _jupyter_labextension_paths():
    return [{"src": data["jupyterlab"].get("outputDir", "static"), "dest": data["name"]}]

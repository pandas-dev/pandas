# Licensed under a 3-clause BSD style license - see LICENSE.rst

import re

import pytest

from asv import util

from . import tools


def test_check(capsys, basic_conf):
    tmpdir, local, conf, machine_file = basic_conf

    # Test check runs (with full benchmark suite)
    with pytest.raises(util.UserError, match="Benchmark suite check failed"):
        tools.run_asv_with_conf(conf, 'check', "--python=same")

    text, err = capsys.readouterr()

    assert re.search(r"params_examples\.track_wrong_number_of_args: call: "
                     r"wrong number of arguments.*: expected 1, has 2", text)
    assert text.count("wrong number of arguments") == 1

import re

from . import tools


def test_profile_python_same(capsys, basic_conf):
    tmpdir, local, conf, machine_file = basic_conf

    # Test Profile can run with python=same
    tools.run_asv_with_conf(conf, 'profile', '--python=same', "time_secondary.track_value",
                            _machine_file=machine_file)
    text, err = capsys.readouterr()

    # time_with_warnings failure case
    assert re.search(r"^\s+1\s+.*time_secondary.*\(track_value\)", text, re.M)

    # Check that it did not clone or install
    assert "Cloning" not in text
    assert "Installing" not in text

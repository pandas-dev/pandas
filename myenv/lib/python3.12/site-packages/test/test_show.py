# Licensed under a 3-clause BSD style license - see LICENSE.rst
import textwrap
from os.path import abspath, dirname, join

from . import tools

BENCHMARK_DIR = abspath(join(dirname(__file__), 'example_results'))
MACHINE_FILE = abspath(join(dirname(__file__), 'asv-machine.json'))


def test_show(capsys, show_fixture):
    conf = show_fixture

    tools.run_asv_with_conf(conf, 'show')
    text, err = capsys.readouterr()
    assert 'py2.7-Cython-numpy1.8' in text
    assert 'py2.7-numpy1.8' in text
    assert 'py2.7-foo-numpy1.8' in text
    assert 'fcf8c079' in text

    tools.run_asv_with_conf(conf, 'show', 'fcf8c079')
    text, err = capsys.readouterr()
    assert "time_ci_small [cheetah/py2.7-numpy1.8]\n  3.00±0s\n\n" in text

    tools.run_asv_with_conf(conf, 'show', 'fcf8c079', '--machine=cheetah',
                            '--bench=time_ci', '--details')
    text, err = capsys.readouterr()
    expected = textwrap.dedent("""
    Commit: fcf8c079

    time_ci_big [cheetah/py2.7-numpy1.8]
      3.00±1s
      ci_99: (1.50s, 3.50s)

    time_ci_small [cheetah/py2.7-numpy1.8]
      3.00±0s
      ci_99: (3.10s, 3.90s)
    """)
    assert text.strip() == expected.strip()


def test_show_durations(capsys, show_fixture):
    conf = show_fixture

    tools.run_asv_with_conf(conf, 'show', '--machine=cheetah', '--durations')
    text, err = capsys.readouterr()
    assert '13dd6571547f8dd87b24c4e29536d33cc4f335c9  6.00s' in text.strip()

    tools.run_asv_with_conf(conf, 'show', '13dd6571', '--machine=cheetah',
                            '--durations')
    text, err = capsys.readouterr()
    expected = textwrap.dedent("""
    Commit: 13dd6571

    Machine    : cheetah
    Environment: py2.7-Cython-numpy1.8

        <setup_cache example:21>  3.00s
        <build>  2.00s
        time_quantity.time_quantity_array_conversion  1.00s

        total duration: 6.00s
    """)
    assert text.strip() == expected.strip()

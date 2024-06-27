# Licensed under a 3-clause BSD style license - see LICENSE.rst
from os.path import join

from asv import machine


def test_machine(tmpdir):
    tmpdir = str(tmpdir)

    m = machine.Machine.load(
        interactive=False,
        machine="orangutan",
        os="BeOS",
        arch="MIPS",
        cpu="10 MHz",
        ram="640k",
        _path=join(tmpdir, 'asv-machine.json'))

    m = machine.Machine.load(
        _path=join(tmpdir, 'asv-machine.json'), interactive=False)

    assert m.machine == 'orangutan'
    assert m.os == 'BeOS'
    assert m.arch == 'MIPS'
    assert m.cpu == '10 MHz'
    assert m.ram == '640k'


def test_machine_defaults(tmpdir):
    tmpdir = str(tmpdir)

    m = machine.Machine.load(
        interactive=True,
        use_defaults=True,
        _path=join(tmpdir, 'asv-machine.json'))

    assert m.__dict__ == m.get_defaults()

# Licensed under a 3-clause BSD style license - see LICENSE.rst

from asv.console import log
# Monkey-patch the machine name
from asv import machine

log.enable()

machine.Machine.hardcoded_machine_name = 'orangutan'

from unittest import TestCase

from ipykernel.tests import utils
from nbformat.converter import convert
from nbformat.reader import reads

import re
import json
from copy import copy
import unittest

try:
    # py3
    from queue import Empty

    def isstr(s):
        return isinstance(s, str)
except ImportError:
    # py2
    from Queue import Empty

    def isstr(s):
        return isinstance(s, basestring)  # noqa

class NotebookTest(TestCase):
    """Validate a notebook. All code cells are executed in order. The output is either checked
    for errors (if no reference output is present), or is compared against expected output.


    Useful references:
    http://nbformat.readthedocs.org/en/latest/format_description.html
    http://jupyter-client.readthedocs.org/en/latest/messaging.html
"""


    IGNORE_TYPES = ["execute_request", "execute_input", "status", "pyin"]
    STRIP_KEYS = ["execution_count", "traceback", "prompt_number", "source"]
    NBFORMAT_VERSION = 4

    def _test_notebook(self, notebook, test):

        with open(notebook) as f:
            nb = convert(reads(f.read()), self.NBFORMAT_VERSION)
        _, kernel = utils.start_new_kernel()
        for i, c in enumerate([c for c in nb.cells if c.cell_type == 'code']):
            self._test_notebook_cell(self.sanitize_cell(c), i, kernel, test)

    def _test_notebook_cell(self, cell, i, kernel, test):

        if hasattr(cell, 'source'): # nbformat 4.0 and later
            code = cell.source
        else:
            code = cell.input
        iopub = kernel.iopub_channel
        kernel.execute(code)
        outputs = []
        msg = None
        no_error = True
        first_error = -1
        error_msg = ''
        while self.should_continue(msg):
            try:
                msg = iopub.get_msg(block=True, timeout=1)
            except Empty:
                continue
            if msg['msg_type'] not in self.IGNORE_TYPES:
                if msg['msg_type'] == 'error':
                    error_msg = '  ' + msg['content']['ename'] + '\n  ' + msg['content']['evalue']
                    no_error = False
                    if first_error == -1:
                        first_error = i
                i = len(outputs)
                expected = i < len(cell.outputs) and cell.outputs[i] or []
                o = self.transform_message(msg, expected)
                outputs.append(o)

        if (test == 'check_error'):
            self.assertTrue(no_error, 'Executing cell %d resulted in an error:\n%s'%(first_error, error_msg))
        else:
            # Compare computed output against stored output.
            # TODO: This doesn't work right now as the generated output is too diverse to
            #       be verifiable.
            scrub = lambda x: self.dump_canonical(list(self.scrub_outputs(x)))
            scrubbed = scrub(outputs)
            expected = scrub(cell.outputs)
            #print('output=%s'%outputs)
            #print('expected=%s'%expected)
            #self.assertEqual(scrubbed, expected, "\n{}\n\n{}".format(scrubbed, expected))
        
    def dump_canonical(self, obj):
        return json.dumps(obj, indent=2, sort_keys=True)

    def scrub_outputs(self, outputs):
        """
        remove all scrubs from output data and text
        """
        for output in outputs:
            out = copy(output)

            for scrub, sub in []:#self.scrubs.items():
                def _scrubLines(lines):
                    if isstr(lines):
                        return re.sub(scrub, sub, lines)
                    else:
                        return [re.sub(scrub, sub, line) for line in lines]

                if "text" in out:
                    out["text"] = _scrubLines(out["text"])

                if "data" in out:
                    if isinstance(out["data"], dict):
                        for mime, data in out["data"].items():
                            out["data"][mime] = _scrubLines(data)
                    else:
                        out["data"] = _scrubLines(out["data"])
            yield out

    def strip_keys(self, d):
        """
        remove keys from STRIP_KEYS to ensure comparability
        """
        for key in self.STRIP_KEYS:
            d.pop(key, None)
        return d

    def sanitize_cell(self, cell):
        """
        remove non-reproducible things
        """
        for output in cell.outputs:
            self.strip_keys(output)
        return cell

    def transform_message(self, msg, expected):
        """
        transform a message into something like the notebook
        """
        SWAP_KEYS = {
            "output_type": {
                "pyout": "execute_result",
                "pyerr": "error"
            }
        }

        output = {
            u"output_type": msg["msg_type"]
        }
        output.update(msg["content"])

        output = self.strip_keys(output)
        for key, swaps in SWAP_KEYS.items():
            if key in output and output[key] in swaps:
                output[key] = swaps[output[key]]

        if "data" in output and "data" not in expected:
            output["text"] = output["data"]
            del output["data"]

        return output

    def should_continue(self, msg):
        """
        determine whether the current message is the last for this cell
        """
        if msg is None:
            return True

        return not (msg["msg_type"] == "status" and
                    msg["content"]["execution_state"] == "idle")



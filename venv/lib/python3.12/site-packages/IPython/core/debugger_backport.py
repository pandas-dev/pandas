"""
The code in this module is a backport of cPython changes in Pdb
that were introduced in Python 3.13 by gh-83151: Make closure work on pdb
https://github.com/python/cpython/pull/111094.
This file should be removed once IPython drops supports for Python 3.12.

The only changes are:
- reformatting by darker (black) formatter
- addition of type-ignore comments to satisfy mypy

Copyright (c) 2001 Python Software Foundation; All Rights Reserved

PYTHON SOFTWARE FOUNDATION LICENSE VERSION 2
--------------------------------------------

1. This LICENSE AGREEMENT is between the Python Software Foundation
("PSF"), and the Individual or Organization ("Licensee") accessing and
otherwise using this software ("Python") in source or binary form and
its associated documentation.

2. Subject to the terms and conditions of this License Agreement, PSF hereby
grants Licensee a nonexclusive, royalty-free, world-wide license to reproduce,
analyze, test, perform and/or display publicly, prepare derivative works,
distribute, and otherwise use Python alone or in any derivative version,
provided, however, that PSF's License Agreement and PSF's notice of copyright,
i.e., "Copyright (c) 2001 Python Software Foundation; All Rights Reserved"
are retained in Python alone or in any derivative version prepared by Licensee.

3. In the event Licensee prepares a derivative work that is based on
or incorporates Python or any part thereof, and wants to make
the derivative work available to others as provided herein, then
Licensee hereby agrees to include in any such work a brief summary of
the changes made to Python.

4. PSF is making Python available to Licensee on an "AS IS"
basis.  PSF MAKES NO REPRESENTATIONS OR WARRANTIES, EXPRESS OR
IMPLIED.  BY WAY OF EXAMPLE, BUT NOT LIMITATION, PSF MAKES NO AND
DISCLAIMS ANY REPRESENTATION OR WARRANTY OF MERCHANTABILITY OR FITNESS
FOR ANY PARTICULAR PURPOSE OR THAT THE USE OF PYTHON WILL NOT
INFRINGE ANY THIRD PARTY RIGHTS.

5. PSF SHALL NOT BE LIABLE TO LICENSEE OR ANY OTHER USERS OF PYTHON
FOR ANY INCIDENTAL, SPECIAL, OR CONSEQUENTIAL DAMAGES OR LOSS AS
A RESULT OF MODIFYING, DISTRIBUTING, OR OTHERWISE USING PYTHON,
OR ANY DERIVATIVE THEREOF, EVEN IF ADVISED OF THE POSSIBILITY THEREOF.

6. This License Agreement will automatically terminate upon a material
breach of its terms and conditions.

7. Nothing in this License Agreement shall be deemed to create any
relationship of agency, partnership, or joint venture between PSF and
Licensee.  This License Agreement does not grant permission to use PSF
trademarks or trade name in a trademark sense to endorse or promote
products or services of Licensee, or any third party.

8. By copying, installing or otherwise using Python, Licensee
agrees to be bound by the terms and conditions of this License
Agreement.
"""

import sys
import types
import codeop
import textwrap
from types import CodeType


class PdbClosureBackport:
    def _exec_in_closure(self, source, globals, locals):  # type: ignore[no-untyped-def]
        """Run source code in closure so code object created within source
        can find variables in locals correctly
        returns True if the source is executed, False otherwise
        """

        # Determine if the source should be executed in closure. Only when the
        # source compiled to multiple code objects, we should use this feature.
        # Otherwise, we can just raise an exception and normal exec will be used.

        code = compile(source, "<string>", "exec")
        if not any(isinstance(const, CodeType) for const in code.co_consts):
            return False

        # locals could be a proxy which does not support pop
        # copy it first to avoid modifying the original locals
        locals_copy = dict(locals)

        locals_copy["__pdb_eval__"] = {"result": None, "write_back": {}}

        # If the source is an expression, we need to print its value
        try:
            compile(source, "<string>", "eval")
        except SyntaxError:
            pass
        else:
            source = "__pdb_eval__['result'] = " + source

        # Add write-back to update the locals
        source = (
            "try:\n"
            + textwrap.indent(source, "  ")
            + "\n"
            + "finally:\n"
            + "  __pdb_eval__['write_back'] = locals()"
        )

        # Build a closure source code with freevars from locals like:
        # def __pdb_outer():
        #   var = None
        #   def __pdb_scope():  # This is the code object we want to execute
        #     nonlocal var
        #     <source>
        #   return __pdb_scope.__code__
        source_with_closure = (
            "def __pdb_outer():\n"
            + "\n".join(f"  {var} = None" for var in locals_copy)
            + "\n"
            + "  def __pdb_scope():\n"
            + "\n".join(f"    nonlocal {var}" for var in locals_copy)
            + "\n"
            + textwrap.indent(source, "    ")
            + "\n"
            + "  return __pdb_scope.__code__"
        )

        # Get the code object of __pdb_scope()
        # The exec fills locals_copy with the __pdb_outer() function and we can call
        # that to get the code object of __pdb_scope()
        ns = {}
        try:
            exec(source_with_closure, {}, ns)
        except Exception:
            return False
        code = ns["__pdb_outer"]()

        cells = tuple(types.CellType(locals_copy.get(var)) for var in code.co_freevars)

        try:
            exec(code, globals, locals_copy, closure=cells)
        except Exception:
            return False

        # get the data we need from the statement
        pdb_eval = locals_copy["__pdb_eval__"]

        # __pdb_eval__ should not be updated back to locals
        pdb_eval["write_back"].pop("__pdb_eval__")

        # Write all local variables back to locals
        locals.update(pdb_eval["write_back"])
        eval_result = pdb_eval["result"]
        if eval_result is not None:
            print(repr(eval_result))

        return True

    def default(self, line):  # type: ignore[no-untyped-def]
        if line[:1] == "!":
            line = line[1:].strip()
        locals = self.curframe_locals
        globals = self.curframe.f_globals
        try:
            buffer = line
            if (
                code := codeop.compile_command(line + "\n", "<stdin>", "single")
            ) is None:
                # Multi-line mode
                with self._disable_command_completion():
                    buffer = line
                    continue_prompt = "...   "
                    while (
                        code := codeop.compile_command(buffer, "<stdin>", "single")
                    ) is None:
                        if self.use_rawinput:
                            try:
                                line = input(continue_prompt)
                            except (EOFError, KeyboardInterrupt):
                                self.lastcmd = ""
                                print("\n")
                                return
                        else:
                            self.stdout.write(continue_prompt)
                            self.stdout.flush()
                            line = self.stdin.readline()
                            if not len(line):
                                self.lastcmd = ""
                                self.stdout.write("\n")
                                self.stdout.flush()
                                return
                            else:
                                line = line.rstrip("\r\n")
                        buffer += "\n" + line
            save_stdout = sys.stdout
            save_stdin = sys.stdin
            save_displayhook = sys.displayhook
            try:
                sys.stdin = self.stdin
                sys.stdout = self.stdout
                sys.displayhook = self.displayhook
                if not self._exec_in_closure(buffer, globals, locals):
                    exec(code, globals, locals)
            finally:
                sys.stdout = save_stdout
                sys.stdin = save_stdin
                sys.displayhook = save_displayhook
        except:
            self._error_exc()

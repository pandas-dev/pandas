"""Miscellaneous inspection tools
"""
from tempfile import NamedTemporaryFile, TemporaryDirectory
import os
import warnings

from numba.core.errors import NumbaWarning


def disassemble_elf_to_cfg(elf, mangled_symbol):
    """
    Gets the CFG of the disassembly of an ELF object, elf, at mangled name,
    mangled_symbol, and renders it appropriately depending on the execution
    environment (terminal/notebook).
    """
    try:
        import r2pipe
    except ImportError:
        raise RuntimeError("r2pipe package needed for disasm CFG")

    def get_rendering(cmd=None):
        from numba.pycc.platform import Toolchain # import local, circular ref
        if cmd is None:
            raise ValueError("No command given")

        with TemporaryDirectory() as tmpdir:
            # Write ELF as a temporary file in the temporary dir, do not delete!
            with NamedTemporaryFile(delete=False, dir=tmpdir) as f:
                f.write(elf)
                f.flush()  # force write, radare2 needs a binary blob on disk

            # Now try and link the ELF, this helps radare2 _a lot_
            linked = False
            try:
                raw_dso_name = f'{os.path.basename(f.name)}.so'
                linked_dso = os.path.join(tmpdir, raw_dso_name)
                tc = Toolchain()
                tc.link_shared(linked_dso, (f.name,))
                obj_to_analyse = linked_dso
                linked = True
            except Exception as e:
                # link failed, mention it to user, radare2 will still be able to
                # analyse the object, but things like dwarf won't appear in the
                # asm as comments.
                msg = ('Linking the ELF object with the distutils toolchain '
                       f'failed with: {e}. Disassembly will still work but '
                       'might be less accurate and will not use DWARF '
                       'information.')
                warnings.warn(NumbaWarning(msg))
                obj_to_analyse = f.name

            # catch if r2pipe can actually talk to radare2
            try:
                flags = ['-2', # close stderr to hide warnings
                         '-e io.cache=true', # fix relocations in disassembly
                         '-e scr.color=1',  # 16 bit ANSI colour terminal
                         '-e asm.dwarf=true', # DWARF decode
                         '-e scr.utf8=true', # UTF8 output looks better
                         ]
                r = r2pipe.open(obj_to_analyse, flags=flags)
                r.cmd('aaaaaa') # analyse as much as possible
                # If the elf is linked then it's necessary to seek as the
                # DSO ctor/dtor is at the default position
                if linked:
                    # r2 only matches up to 61 chars?! found this by experiment!
                    mangled_symbol_61char = mangled_symbol[:61]
                    # switch off demangle, the seek is on a mangled symbol
                    r.cmd('e bin.demangle=false')
                    # seek to the mangled symbol address
                    r.cmd(f's `is~ {mangled_symbol_61char}[1]`')
                    # switch demangling back on for output purposes
                    r.cmd('e bin.demangle=true')
                data = r.cmd('%s' % cmd) # print graph
                r.quit()
            except Exception as e:
                if "radare2 in PATH" in str(e):
                    msg = ("This feature requires 'radare2' to be "
                           "installed and available on the system see: "
                           "https://github.com/radareorg/radare2. "
                           "Cannot find 'radare2' in $PATH.")
                    raise RuntimeError(msg)
                else:
                    raise e
        return data

    class DisasmCFG(object):

        def _repr_svg_(self):
            try:
                import graphviz
            except ImportError:
                raise RuntimeError("graphviz package needed for disasm CFG")
            jupyter_rendering = get_rendering(cmd='agfd')
            # this just makes it read slightly better in jupyter notebooks
            jupyter_rendering.replace('fontname="Courier",',
                                      'fontname="Courier",fontsize=6,')
            src = graphviz.Source(jupyter_rendering)
            return src.pipe('svg').decode('UTF-8')

        def __repr__(self):
            return get_rendering(cmd='agf')

    return DisasmCFG()

from llvmlite.ir import CallInstr


class Visitor(object):
    def visit(self, module):
        self._module = module
        for func in module.functions:
            self.visit_Function(func)

    def visit_Function(self, func):
        self._function = func
        for bb in func.blocks:
            self.visit_BasicBlock(bb)

    def visit_BasicBlock(self, bb):
        self._basic_block = bb
        for instr in bb.instructions:
            self.visit_Instruction(instr)

    def visit_Instruction(self, instr):
        raise NotImplementedError

    @property
    def module(self):
        return self._module

    @property
    def function(self):
        return self._function

    @property
    def basic_block(self):
        return self._basic_block


class CallVisitor(Visitor):
    def visit_Instruction(self, instr):
        if isinstance(instr, CallInstr):
            self.visit_Call(instr)

    def visit_Call(self, instr):
        raise NotImplementedError


class ReplaceCalls(CallVisitor):
    def __init__(self, orig, repl):
        super(ReplaceCalls, self).__init__()
        self.orig = orig
        self.repl = repl
        self.calls = []

    def visit_Call(self, instr):
        if instr.callee == self.orig:
            instr.replace_callee(self.repl)
            self.calls.append(instr)


def replace_all_calls(mod, orig, repl):
    """Replace all calls to `orig` to `repl` in module `mod`.
    Returns the references to the returned calls
    """
    rc = ReplaceCalls(orig, repl)
    rc.visit(mod)
    return rc.calls

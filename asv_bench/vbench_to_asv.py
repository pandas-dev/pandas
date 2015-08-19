import ast
import vbench
import os
import sys
import astor
import glob


def vbench_to_asv_source(bench, kinds=None):
    tab = ' ' * 4
    if kinds is None:
        kinds = ['time']

    output = 'class {}(object):\n'.format(bench.name)
    output += tab + 'goal_time = 0.2\n\n'

    if bench.setup:
        indented_setup = [tab * 2 + '{}\n'.format(x) for x in bench.setup.splitlines()]
        output += tab + 'def setup(self):\n' + ''.join(indented_setup) + '\n'

    for kind in kinds:
        output += tab + 'def {}_{}(self):\n'.format(kind, bench.name)
        for line in bench.code.splitlines():
            output += tab * 2 + line + '\n'
        output += '\n\n'

    if bench.cleanup:
        output += tab + 'def teardown(self):\n' + tab * 2 + bench.cleanup

    output += '\n\n'
    return output


class AssignToSelf(ast.NodeTransformer):
    def __init__(self):
        super(AssignToSelf, self).__init__()
        self.transforms = {}
        self.imports = []

        self.in_class_define = False
        self.in_setup = False

    def visit_ClassDef(self, node):
        self.transforms = {}
        self.in_class_define = True
        self.generic_visit(node)
        return node

    def visit_TryExcept(self, node):
        if any([isinstance(x, (ast.Import, ast.ImportFrom)) for x in node.body]):
            self.imports.append(node)
        else:
            self.generic_visit(node)
            return node

    def visit_Assign(self, node):
        for target in node.targets:
            if isinstance(target, ast.Name) and not isinstance(target.ctx, ast.Param) and not self.in_class_define:
                self.transforms[target.id] = 'self.' + target.id
        self.generic_visit(node)

        return node

    def visit_Name(self, node):
        new_node = node
        if node.id in self.transforms:
            if not isinstance(node.ctx, ast.Param):
                new_node = ast.Attribute(value=ast.Name(id='self', ctx=node.ctx), attr=node.id, ctx=node.ctx)

        self.generic_visit(node)

        return ast.copy_location(new_node, node)

    def visit_Import(self, node):
        self.imports.append(node)

    def visit_ImportFrom(self, node):
        self.imports.append(node)

    def visit_FunctionDef(self, node):
        """Delete functions that are empty due to imports being moved"""
        self.in_class_define = False

        if self.in_setup:
            node.col_offset -= 4
            ast.increment_lineno(node, -1)

        if node.name == 'setup':
            self.in_setup = True

        self.generic_visit(node)

        if node.name == 'setup':
            self.in_setup = False

        if node.body:
            return node


def translate_module(target_module):
    g_vars = {}
    l_vars = {}
    exec('import ' + target_module) in g_vars

    print target_module
    module = eval(target_module, g_vars)

    benchmarks = []
    for obj_str in dir(module):
        obj = getattr(module, obj_str)
        if isinstance(obj, vbench.benchmark.Benchmark):
            benchmarks.append(obj)

    if not benchmarks:
        return

    rewritten_output = ''
    for bench in benchmarks:
        rewritten_output += vbench_to_asv_source(bench)

    with open('rewrite.py', 'w') as f:
        f.write(rewritten_output)

    ast_module = ast.parse(rewritten_output)

    transformer = AssignToSelf()
    transformed_module = transformer.visit(ast_module)

    unique_imports = {astor.to_source(node): node for node in transformer.imports}

    transformed_module.body = unique_imports.values() + transformed_module.body

    transformed_source = astor.to_source(transformed_module)

    with open('benchmarks/{}.py'.format(target_module), 'w') as f:
        f.write(transformed_source)


if __name__ == '__main__':
    cwd = os.getcwd()
    new_dir = os.path.join(os.path.dirname(__file__), '../vb_suite')
    sys.path.insert(0, new_dir)

    for module in glob.glob(os.path.join(new_dir, '*.py')):
        mod = os.path.basename(module)
        if mod in ['make.py', 'measure_memory_consumption.py', 'perf_HEAD.py', 'run_suite.py', 'test_perf.py', 'generate_rst_files.py', 'test.py', 'suite.py']:
            continue
        print
        print mod

        translate_module(mod.replace('.py', ''))

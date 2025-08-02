#######################################################################
# Copyright (c) 2019-present, Blosc Development Team <blosc@blosc.org>
# All rights reserved.
#
# This source code is licensed under a BSD-style license (found in the
# LICENSE file in the root directory of this source tree)
#######################################################################

import re


def _inherit_doc_parameter(parent_func, parameter, replacements=None):
    # Small decorator to copy parameter descriptions from other functions with optional support of replacement rules
    # (see blosc2.open for an example)
    def wrapper(child_func):
        # Copy relevant lines from parent function
        matching_lines = []
        indent_parent = None
        for line in parent_func.__doc__.splitlines():
            if parameter in line:
                match = re.search(rf"(\s*){parameter}", line)
                assert match is not None, (
                    f"Parameter {parameter} not found in the docstring of {parent_func.__name__}"
                )
                indent_parent = match.group(1)

                # The first line should be without the indentation because it will be placed at the correct location
                # in the child function
                matching_lines.append(line.lstrip())
            elif indent_parent is not None:
                if re.search(rf"^{indent_parent}\w+", line) is not None:
                    # Next parameter starts, stop copying lines
                    break
                matching_lines.append(line)
        assert matching_lines, (
            f"Could not extract the parameter {parameter} from the docstring of {parent_func.__name__}"
        )

        # Replace the indentation of the parent with the indentation used in the child function
        match = re.search(rf"([ \t]+){parameter}", child_func.__doc__)
        assert match is not None, (
            f"Parameter {parameter} not found in the docstring of {child_func.__name__}"
        )
        indent_child = match.group(1)

        # First line contains the parameter name itself which should not be indented
        matching_lines = [matching_lines[0].lstrip()] + [
            ml.replace(indent_parent, indent_child, 1) for ml in matching_lines[1:]
        ]

        child_func.__doc__ = child_func.__doc__.replace(parameter, "\n".join(matching_lines))

        if replacements is not None:
            for regex, repl in replacements.items():
                new_doc = re.sub(regex, repl, child_func.__doc__)
                assert new_doc != child_func.__doc__, (
                    f"The replacement rule {regex}: {repl} did not change the docstring of {child_func.__name__}"
                )
                child_func.__doc__ = new_doc

        return child_func

    return wrapper

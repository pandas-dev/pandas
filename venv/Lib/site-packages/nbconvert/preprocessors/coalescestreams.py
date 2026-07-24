"""Preprocessor for merging consecutive stream outputs for easier handling."""

import re

# Copyright (c) IPython Development Team.
# Distributed under the terms of the Modified BSD License.
from nbconvert.preprocessors import Preprocessor

CR_PAT = re.compile(r".*\r(?=[^\n])")


class CoalesceStreamsPreprocessor(Preprocessor):
    """
    Merge consecutive sequences of stream output into single stream
    to prevent extra newlines inserted at flush calls
    """

    def preprocess_cell(self, cell, resources, cell_index):
        """
        Apply a transformation on each cell. See base.py for details.
        """
        outputs = cell.get("outputs", [])
        if not outputs:
            return cell, resources

        last = outputs[0]
        new_outputs = [last]
        for output in outputs[1:]:
            if (
                output.output_type == "stream"
                and last.output_type == "stream"
                and last.name == output.name
            ):
                last.text += output.text
            else:
                new_outputs.append(output)
                last = output

        # process \r characters
        for output in new_outputs:
            if output.output_type == "stream" and "\r" in output.text:
                output.text = CR_PAT.sub("", output.text)

        cell.outputs = new_outputs
        return cell, resources

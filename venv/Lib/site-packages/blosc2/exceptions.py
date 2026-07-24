#######################################################################
# Copyright (c) 2019-present, Blosc Development Team <blosc@blosc.org>
# All rights reserved.
#
# SPDX-License-Identifier: BSD-3-Clause
#######################################################################


class MissingOperands(ValueError):
    def __init__(self, expr, missing_ops):
        self.expr = expr
        self.missing_ops = missing_ops

        message = f'Lazy expression "{expr}" with missing operands: {missing_ops}'
        super().__init__(message)

from __future__ import annotations

import math

import numpy as np

from dask.array._shuffle import _calculate_new_chunksizes
from dask.array.core import asarray, blockwise, einsum_lookup
from dask.utils import cached_max, derived_from

einsum_symbols = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ"
einsum_symbols_set = set(einsum_symbols)
from dask import config


def chunk_einsum(*operands, **kwargs):
    subscripts = kwargs.pop("subscripts")
    ncontract_inds = kwargs.pop("ncontract_inds")
    dtype = kwargs.pop("kernel_dtype")
    einsum = einsum_lookup.dispatch(type(operands[0]))
    chunk = einsum(subscripts, *operands, dtype=dtype, **kwargs)

    # Avoid concatenate=True in blockwise by adding 1's
    # for the contracted dimensions
    return chunk.reshape(chunk.shape + (1,) * ncontract_inds)


# This function duplicates numpy's _parse_einsum_input() function
# See https://github.com/numpy/numpy/blob/master/LICENSE.txt
# or NUMPY_LICENSE.txt within this directory
def parse_einsum_input(operands):
    """
    A reproduction of numpy's _parse_einsum_input()
    which in itself is a reproduction of
    c side einsum parsing in python.

    Returns
    -------
    input_strings : str
        Parsed input strings
    output_string : str
        Parsed output string
    operands : list of array_like
        The operands to use in the numpy contraction
    Examples
    --------
    The operand list is simplified to reduce printing:
    >> a = np.random.rand(4, 4)
    >> b = np.random.rand(4, 4, 4)
    >> __parse_einsum_input(('...a,...a->...', a, b))
    ('za,xza', 'xz', [a, b])
    >> __parse_einsum_input((a, [Ellipsis, 0], b, [Ellipsis, 0]))
    ('za,xza', 'xz', [a, b])
    """

    if len(operands) == 0:
        raise ValueError("No input operands")

    if isinstance(operands[0], str):
        subscripts = operands[0].replace(" ", "")
        operands = [asarray(o) for o in operands[1:]]

        # Ensure all characters are valid
        for s in subscripts:
            if s in ".,->":
                continue
            if s not in einsum_symbols_set:
                raise ValueError("Character %s is not a valid symbol." % s)

    else:
        tmp_operands = list(operands)
        operand_list = []
        subscript_list = []
        for _ in range(len(operands) // 2):
            operand_list.append(tmp_operands.pop(0))
            subscript_list.append(tmp_operands.pop(0))

        output_list = tmp_operands[-1] if len(tmp_operands) else None
        operands = [asarray(v) for v in operand_list]
        subscripts = ""
        last = len(subscript_list) - 1
        for num, sub in enumerate(subscript_list):
            for s in sub:
                if s is Ellipsis:
                    subscripts += "..."
                elif isinstance(s, int):
                    subscripts += einsum_symbols[s]
                else:
                    raise TypeError(
                        "For this input type lists must contain "
                        "either int or Ellipsis"
                    )
            if num != last:
                subscripts += ","

        if output_list is not None:
            subscripts += "->"
            for s in output_list:
                if s is Ellipsis:
                    subscripts += "..."
                elif isinstance(s, int):
                    subscripts += einsum_symbols[s]
                else:
                    raise TypeError(
                        "For this input type lists must contain "
                        "either int or Ellipsis"
                    )
    # Check for proper "->"
    if ("-" in subscripts) or (">" in subscripts):
        invalid = (subscripts.count("-") > 1) or (subscripts.count(">") > 1)
        if invalid or (subscripts.count("->") != 1):
            raise ValueError("Subscripts can only contain one '->'.")

    # Parse ellipses
    if "." in subscripts:
        used = subscripts.replace(".", "").replace(",", "").replace("->", "")
        unused = list(einsum_symbols_set - set(used))
        ellipse_inds = "".join(unused)
        longest = 0

        if "->" in subscripts:
            input_tmp, output_sub = subscripts.split("->")
            split_subscripts = input_tmp.split(",")
            out_sub = True
        else:
            split_subscripts = subscripts.split(",")
            out_sub = False

        for num, sub in enumerate(split_subscripts):
            if "." in sub:
                if (sub.count(".") != 3) or (sub.count("...") != 1):
                    raise ValueError("Invalid Ellipses.")

                # Take into account numerical values
                if operands[num].shape == ():
                    ellipse_count = 0
                else:
                    ellipse_count = max(operands[num].ndim, 1)
                    ellipse_count -= len(sub) - 3

                if ellipse_count > longest:
                    longest = ellipse_count

                if ellipse_count < 0:
                    raise ValueError("Ellipses lengths do not match.")
                elif ellipse_count == 0:
                    split_subscripts[num] = sub.replace("...", "")
                else:
                    rep_inds = ellipse_inds[-ellipse_count:]
                    split_subscripts[num] = sub.replace("...", rep_inds)

        subscripts = ",".join(split_subscripts)
        if longest == 0:
            out_ellipse = ""
        else:
            out_ellipse = ellipse_inds[-longest:]

        if out_sub:
            subscripts += "->" + output_sub.replace("...", out_ellipse)
        else:
            # Special care for outputless ellipses
            output_subscript = ""
            tmp_subscripts = subscripts.replace(",", "")
            for s in sorted(set(tmp_subscripts)):
                if s not in einsum_symbols_set:
                    raise ValueError("Character %s is not a valid symbol." % s)
                if tmp_subscripts.count(s) == 1:
                    output_subscript += s
            normal_inds = "".join(sorted(set(output_subscript) - set(out_ellipse)))

            subscripts += "->" + out_ellipse + normal_inds

    # Build output string if does not exist
    if "->" in subscripts:
        input_subscripts, output_subscript = subscripts.split("->")
    else:
        input_subscripts = subscripts
        # Build output subscripts
        tmp_subscripts = subscripts.replace(",", "")
        output_subscript = ""
        for s in sorted(set(tmp_subscripts)):
            if s not in einsum_symbols_set:
                raise ValueError("Character %s is not a valid symbol." % s)
            if tmp_subscripts.count(s) == 1:
                output_subscript += s

    # Make sure output subscripts are in the input
    for char in output_subscript:
        if char not in input_subscripts:
            raise ValueError("Output character %s did not appear in the input" % char)

    # Make sure number operands is equivalent to the number of terms
    if len(input_subscripts.split(",")) != len(operands):
        raise ValueError(
            "Number of einsum subscripts must be equal to the number of operands."
        )

    return (input_subscripts, output_subscript, operands)


@derived_from(np)
def einsum(*operands, dtype=None, optimize=False, split_every=None, **kwargs):
    """Dask added an additional keyword-only argument ``split_every``.

    split_every: int >= 2 or dict(axis: int), optional
        Determines the depth of the recursive aggregation.
        Defaults to ``None`` which would let dask heuristically
        decide a good default.
    """

    einsum_dtype = dtype

    inputs, outputs, ops = parse_einsum_input(operands)
    subscripts = "->".join((inputs, outputs))

    # Infer the output dtype from operands
    if dtype is None:
        dtype = np.result_type(*[o.dtype for o in ops])

    if optimize is not False:
        # Avoid computation of dask arrays within np.einsum_path
        # by passing in small numpy arrays broadcasted
        # up to the right shape
        fake_ops = [np.broadcast_to(o.dtype.type(0), shape=o.shape) for o in ops]
        optimize, _ = np.einsum_path(subscripts, *fake_ops, optimize=optimize)

    inputs = [tuple(i) for i in inputs.split(",")]

    # Set of all indices
    all_inds = {a for i in inputs for a in i}

    # Which indices are contracted?
    contract_inds = all_inds - set(outputs)
    ncontract_inds = len(contract_inds)

    if len(inputs) > 1 and len(outputs) > 0:
        # Calculate the increase in chunk size compared to the largest input chunk
        max_chunk_sizes, max_chunk_size_input = {}, 1
        for op, input in zip(ops, inputs):
            max_chunk_size_input = max(
                math.prod(map(cached_max, op.chunks)), max_chunk_size_input
            )
            max_chunk_sizes.update(
                {
                    inp: max(cached_max(op.chunks[i]), max_chunk_sizes.get(inp, 1))
                    for i, inp in enumerate(input)
                    if inp not in contract_inds
                }
            )

        max_chunk_size_output = math.prod(max_chunk_sizes.values())
        factor = max_chunk_size_output / (
            max_chunk_size_input * config.get("array.chunk-size-tolerance")
        )

        # Rechunk inputs to make input chunks smaller to avoid an increase in
        # output chunks
        new_ops = []
        for op, input in zip(ops, inputs):
            changeable_dimensions = {ctr for ctr, i in enumerate(input) if i in outputs}
            f = max(factor ** (len(changeable_dimensions) / len(outputs)), 1)
            result = _calculate_new_chunksizes(
                op.chunks,
                list(op.chunks),
                changeable_dimensions,
                math.prod(map(cached_max, op.chunks)) / f,
            )
            new_ops.append(op.rechunk(result))
        ops = new_ops

    # Introduce the contracted indices into the blockwise product
    # so that we get numpy arrays, not lists
    result = blockwise(
        chunk_einsum,
        tuple(outputs) + tuple(contract_inds),
        *(a for ap in zip(ops, inputs) for a in ap),
        # blockwise parameters
        adjust_chunks={ind: 1 for ind in contract_inds},
        dtype=dtype,
        # np.einsum parameters
        subscripts=subscripts,
        kernel_dtype=einsum_dtype,
        ncontract_inds=ncontract_inds,
        optimize=optimize,
        **kwargs,
    )

    # Now reduce over any extra contraction dimensions
    if ncontract_inds > 0:
        size = len(outputs)
        return result.sum(
            axis=list(range(size, size + ncontract_inds)), split_every=split_every
        )

    return result

from io import StringIO

import numpy as np
import pytest

import pandas as pd
import pandas._testing as tm


class BaseParsingTests:
    @pytest.mark.parametrize("engine", ["c", "python"])
    def test_EA_types(self, engine, data, request):
        if engine == "c" and data.dtype.kind == "c":
            request.node.add_marker(
                pytest.mark.xfail(
                    reason=f"engine '{engine}' cannot parse the dtype {data.dtype.name}"
                )
            )
        df = pd.DataFrame({"with_dtype": pd.Series(data, dtype=str(data.dtype))})
        csv_output = df.to_csv(index=False, na_rep=np.nan)
        result = pd.read_csv(
            StringIO(csv_output), dtype={"with_dtype": str(data.dtype)}, engine=engine
        )
        expected = df
        tm.assert_frame_equal(result, expected)

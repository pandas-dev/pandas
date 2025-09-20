from numba import cuda
from numba.cuda.testing import CUDATestCase, skip_on_cudasim, skip_unless_cc_53


class TestIsFP16Supported(CUDATestCase):
    def test_is_fp16_supported(self):
        self.assertTrue(cuda.is_float16_supported())

    @skip_on_cudasim
    @skip_unless_cc_53
    def test_device_supports_float16(self):
        self.assertTrue(cuda.get_current_device().supports_float16)

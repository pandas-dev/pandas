from numba.cuda.deviceufunc import GUFuncEngine
import unittest


def template(signature, shapes, expects):
    gufb = GUFuncEngine.from_signature(signature)
    sch = gufb.schedule(shapes)
    for k, v in expects.items():
        got = getattr(sch, k)
        if got != v:
            fmt = 'error for %s: got=%s but expect=%s'
            raise AssertionError(fmt % (k, got, v))


class TestGUFuncScheduling(unittest.TestCase):
    def test_signature_1(self):
        signature = '(m, n), (n, p) -> (m, p)'
        shapes = (100, 4, 5), (1, 5, 7)
        expects = dict(
            ishapes=[(4, 5), (5, 7)],
            oshapes=[(4, 7)],
            loopdims=(100,),
            pinned=[False, True]
        )
        template(signature, shapes, expects)

    def test_signature_2(self):
        signature = '(m, n), (n, p) -> (m, p)'
        shapes = (100, 4, 5), (100, 5, 7)
        expects = dict(
            ishapes=[(4, 5), (5, 7)],
            oshapes=[(4, 7)],
            loopdims=(100,),
            pinned=[False, False]
        )
        template(signature, shapes, expects)

    def test_signature_3(self):
        signature = '(m, n), (n, p) -> (m, p)'
        shapes = (12, 34, 4, 5), (12, 34, 5, 7)
        expects = dict(
            ishapes=[(4, 5), (5, 7)],
            oshapes=[(4, 7)],
            loopdims=(12, 34),
            pinned=[False, False]
        )
        template(signature, shapes, expects)

    def test_signature_4(self):
        signature = '(m, n), (n, p) -> (m, p)'
        shapes = (4, 5), (5, 7)
        expects = dict(
            ishapes=[(4, 5), (5, 7)],
            oshapes=[(4, 7)],
            loopdims=(),
            pinned=[False, False]
        )
        template(signature, shapes, expects)

    def test_signature_5(self):
        signature = '(a), (a) -> (a)'
        shapes = (5,), (5,)
        expects = dict(
            ishapes=[(5,), (5,)],
            oshapes=[(5,)],
            loopdims=(),
            pinned=[False, False]
        )
        template(signature, shapes, expects)

    def test_signature_6(self):
        signature = '(), () -> ()'
        shapes = (5,), (5,)
        expects = dict(
            ishapes=[(), ()],
            oshapes=[()],
            loopdims=(5,),
            pinned=[False, False]
        )
        template(signature, shapes, expects)

    def test_signature_7(self):
        signature = '(), () -> ()'
        shapes = (5,), ()
        expects = dict(
            ishapes=[(), ()],
            oshapes=[()],
            loopdims=(5,),
            pinned=[False, True]
        )
        template(signature, shapes, expects)


if __name__ == '__main__':
    unittest.main()

# Simple support library for our run tests.

from contextlib import contextmanager
from typing import Iterator, TypeVar, Generator, Optional, List, Tuple, Sequence, Union

@contextmanager
def assertRaises(typ: type, msg: str = '') -> Iterator[None]:
    try:
        yield
    except Exception as e:
        assert isinstance(e, typ), "{} is not a {}".format(e, typ.__name__)
        assert msg in str(e), 'Message "{}" does not match "{}"'.format(e, msg)
    else:
        assert False, "Expected {} but got no exception".format(typ.__name__)

T = TypeVar('T')
U = TypeVar('U')
V = TypeVar('V')

def run_generator(gen: Generator[T, V, U],
                  inputs: Optional[List[V]] = None,
                  p: bool = False) -> Tuple[Sequence[T], Union[U, str]]:
    res = []  # type: List[T]
    i = -1
    while True:
        try:
            if i >= 0 and inputs:
                # ... fixtures don't have send
                val = gen.send(inputs[i])  # type: ignore
            else:
                val = next(gen)
        except StopIteration as e:
            return (tuple(res), e.value)
        except Exception as e:
            return (tuple(res), str(e))
        if p:
            print(val)
        res.append(val)
        i += 1

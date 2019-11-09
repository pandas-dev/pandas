from dataclasses import dataclass
# cd vortexa_utils/
# from aws.utils.dataclasses import nested_dataclass
from vortexa_utils.aws.utils.dataclasses import nested_dataclass


@dataclass
class Foo:
    a: str
    b: int


@nested_dataclass
class Bar:
    foo: Foo
    baz: str


@nested_dataclass
class Bill:
    bar: Bar


def test_init_class():
    data = dict(
        bar=dict(
            foo=dict(a="hello", b=1),
            baz="world"
        )
    )
    foo = Foo(**data['bar']['foo'])
    bar = Bar(**data['bar'])
    bill = Bill(**data)

    assert bill.bar == bar
    assert bill.bar.foo == foo

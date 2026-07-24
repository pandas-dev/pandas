from _typeshed import Incomplete

from pony.orm.core import Database, Entity

db: Database

class Customer(Entity):
    __slots__ = ()
    id: Incomplete
    name: Incomplete
    email: Incomplete
    orders: Incomplete

class Order(Entity):
    __slots__ = ()
    id: Incomplete
    total_price: Incomplete
    customer: Incomplete
    items: Incomplete

class Product(Entity):
    __slots__ = ()
    id: Incomplete
    name: Incomplete
    price: Incomplete
    items: Incomplete

class OrderItem(Entity):
    __slots__ = ()
    quantity: Incomplete
    order: Incomplete
    product: Incomplete

def populate_database() -> None: ...

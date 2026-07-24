from _typeshed import Incomplete

from pony.orm.core import Database, Entity

db: Database

class Customer(Entity):
    __slots__ = ()
    email: Incomplete
    password: Incomplete
    name: Incomplete
    country: Incomplete
    address: Incomplete
    cart_items: Incomplete
    orders: Incomplete

class Product(Entity):
    __slots__ = ()
    id: Incomplete
    name: Incomplete
    categories: Incomplete
    description: Incomplete
    picture: Incomplete
    price: Incomplete
    quantity: Incomplete
    cart_items: Incomplete
    order_items: Incomplete

class CartItem(Entity):
    __slots__ = ()
    quantity: Incomplete
    customer: Incomplete
    product: Incomplete

class OrderItem(Entity):
    __slots__ = ()
    quantity: Incomplete
    price: Incomplete
    order: Incomplete
    product: Incomplete

class Order(Entity):
    __slots__ = ()
    id: Incomplete
    state: Incomplete
    date_created: Incomplete
    date_shipped: Incomplete
    date_delivered: Incomplete
    total_price: Incomplete
    customer: Incomplete
    items: Incomplete

class Category(Entity):
    __slots__ = ()
    name: Incomplete
    products: Incomplete

CREATED: str
SHIPPED: str
DELIVERED: str
CANCELLED: str

def populate_database() -> None: ...
def test_queries() -> None: ...

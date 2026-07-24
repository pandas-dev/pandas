from _typeshed import Incomplete

from pony.orm.core import Database, Entity

db: Database

class Person(Entity):
    name: str
    age: int

p1: Person
p2: Person
x: int
y: int
q: Incomplete

from datetime import datetime

from pony.orm.core import Database, Entity

db: Database

class User(Entity):
    login: str
    password: str
    last_login: datetime | None

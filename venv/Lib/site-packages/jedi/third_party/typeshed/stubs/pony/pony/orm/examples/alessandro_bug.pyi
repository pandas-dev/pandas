from _typeshed import Incomplete

from pony.orm import *
from pony.orm.core import Database, Entity

database: Database

class User(Entity):
    __slots__ = ()
    user_id: PrimaryKey
    owned_pokemons: Incomplete
    is_admin: Incomplete
    is_banned: Incomplete
    @staticmethod
    def get_or_create(user_id: int) -> User: ...
    @staticmethod
    def get_by_id(user_id: int) -> User: ...
    def catch_pokemon(self, pokemon: Pokemon): ...
    def remove_pokemon(self, pokemon: Pokemon): ...
    favorite_color: Incomplete
    def set_favorite_color(self, color: tuple[Incomplete, ...]): ...

class Pokemon(Entity):
    __slots__ = ()
    name: Incomplete
    pokemon_id: Incomplete
    sprite: Incomplete
    is_shiny: Incomplete
    owner: Incomplete
    spawned_chat_id: Incomplete
    spawned_message_id: Incomplete
    @property
    def captured(self) -> bool: ...
    def caught_by(self, user: User): ...

class Chat(Entity):
    __slots__ = ()
    chat_id: Incomplete
    active: Incomplete
    def activate(self) -> None: ...
    def deactivate(self) -> None: ...
    @staticmethod
    def get_or_create(chat_id: int) -> Chat: ...
    @staticmethod
    def get_by_id(chat_id: int) -> Chat: ...

def spawn_pokemon(
    chat_id: int, message_id: int, pokemon_json: dict[Incomplete, Incomplete], is_shiny: bool = False
) -> Pokemon: ...
def get_spawned_pokemon(chat_id: int, message_id: int) -> Pokemon | None: ...
def setup() -> None: ...

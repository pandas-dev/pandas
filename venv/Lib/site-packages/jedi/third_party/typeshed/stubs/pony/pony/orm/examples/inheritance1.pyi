from _typeshed import Incomplete

from pony.orm.core import Database, Entity

db: Database

class Person(Entity):
    __slots__ = ()
    id: Incomplete
    name: Incomplete
    dob: Incomplete
    ssn: Incomplete

class Student(Person):
    __slots__ = ()
    group: Incomplete
    mentor: Incomplete
    attend_courses: Incomplete

class Teacher(Person):
    __slots__ = ()
    teach_courses: Incomplete
    apprentices: Incomplete
    salary: Incomplete

class Assistant(Student, Teacher):
    __slots__ = ()

class Professor(Teacher):
    __slots__ = ()
    position: Incomplete

class Group(Entity):
    __slots__ = ()
    number: Incomplete
    students: Incomplete

class Course(Entity):
    __slots__ = ()
    name: Incomplete
    semester: Incomplete
    students: Incomplete
    teachers: Incomplete

def populate_database() -> None: ...
def show_all_persons() -> None: ...

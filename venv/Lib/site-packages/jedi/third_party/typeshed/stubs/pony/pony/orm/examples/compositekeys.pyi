from _typeshed import Incomplete

from pony.orm.core import Database, Entity

db: Database

class Group(Entity):
    __slots__ = ()
    dept: Incomplete
    year: Incomplete
    spec: Incomplete
    students: Incomplete
    courses: Incomplete
    lessons: Incomplete

class Department(Entity):
    __slots__ = ()
    number: Incomplete
    faculty: Incomplete
    name: Incomplete
    groups: Incomplete
    teachers: Incomplete

class Faculty(Entity):
    __slots__ = ()
    number: Incomplete
    name: Incomplete
    depts: Incomplete

class Student(Entity):
    __slots__ = ()
    name: Incomplete
    group: Incomplete
    dob: Incomplete
    grades: Incomplete

class Grade(Entity):
    __slots__ = ()
    student: Incomplete
    task: Incomplete
    date: Incomplete
    value: Incomplete

class Task(Entity):
    __slots__ = ()
    course: Incomplete
    type: Incomplete
    number: Incomplete
    descr: Incomplete
    grades: Incomplete

class Course(Entity):
    __slots__ = ()
    subject: Incomplete
    semester: Incomplete
    groups: Incomplete
    tasks: Incomplete
    lessons: Incomplete
    teachers: Incomplete

class Subject(Entity):
    __slots__ = ()
    name: Incomplete
    descr: Incomplete
    courses: Incomplete

class Room(Entity):
    __slots__ = ()
    building: Incomplete
    number: Incomplete
    floor: Incomplete
    schedules: Incomplete

class Teacher(Entity):
    __slots__ = ()
    dept: Incomplete
    name: Incomplete
    courses: Incomplete
    lessons: Incomplete

class Lesson(Entity):
    __slots__ = ()
    groups: Incomplete
    course: Incomplete
    room: Incomplete
    teacher: Incomplete
    date: Incomplete

def test_queries() -> None: ...

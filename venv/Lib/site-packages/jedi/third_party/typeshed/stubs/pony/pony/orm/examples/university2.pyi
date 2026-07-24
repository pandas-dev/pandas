from _typeshed import Incomplete

from pony.orm.core import Database, Entity

db: Database

class Faculty(Entity):
    __slots__ = ()
    number: Incomplete
    name: Incomplete
    departments: Incomplete

class Department(Entity):
    __slots__ = ()
    number: Incomplete
    name: Incomplete
    faculty: Incomplete
    teachers: Incomplete
    majors: Incomplete
    groups: Incomplete

class Group(Entity):
    __slots__ = ()
    number: Incomplete
    grad_year: Incomplete
    department: Incomplete
    lessons: Incomplete
    students: Incomplete

class Student(Entity):
    __slots__ = ()
    name: Incomplete
    scholarship: Incomplete
    group: Incomplete
    grades: Incomplete

class Major(Entity):
    __slots__ = ()
    name: Incomplete
    department: Incomplete
    courses: Incomplete

class Subject(Entity):
    __slots__ = ()
    name: Incomplete
    courses: Incomplete
    teachers: Incomplete

class Course(Entity):
    __slots__ = ()
    major: Incomplete
    subject: Incomplete
    semester: Incomplete
    lect_hours: Incomplete
    pract_hours: Incomplete
    credit: Incomplete
    lessons: Incomplete
    grades: Incomplete

class Lesson(Entity):
    __slots__ = ()
    day_of_week: Incomplete
    meeting_time: Incomplete
    classroom: Incomplete
    course: Incomplete
    teacher: Incomplete
    groups: Incomplete

class Grade(Entity):
    __slots__ = ()
    student: Incomplete
    course: Incomplete
    teacher: Incomplete
    date: Incomplete
    value: Incomplete

class Teacher(Entity):
    __slots__ = ()
    name: Incomplete
    degree: Incomplete
    department: Incomplete
    subjects: Incomplete
    lessons: Incomplete
    grades: Incomplete

class Building(Entity):
    __slots__ = ()
    number: Incomplete
    description: Incomplete
    classrooms: Incomplete

class Classroom(Entity):
    __slots__ = ()
    building: Incomplete
    number: Incomplete
    description: Incomplete
    lessons: Incomplete

def test_queries() -> None: ...

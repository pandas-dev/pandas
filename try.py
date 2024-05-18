def process_user(name: str, age: int) -> str:
    return f"{name} is {age} years old."


print(
    process_user("John", "twenty")
)  # mypy will not catch this type error, but Python will run it without complaints

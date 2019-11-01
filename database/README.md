# Vortexa Utils DatabaseFactory

Small factory class to give you a `SqlAlchemy` engine connection to an
`AWS rds` instance ensuring SSL and credentials are obtained with the secrets manager
## Usage

```python
db_factory = DatabaseFactory()
engine = db_factory.engine(dbname='rolling_backup')

sql = """
SELECT
    name
FROM new_polygons where name is not Null;
"""

engine.execute(sql)
```
## TODO Other utility functions

- [ ]  create a `~/.dbpass` file

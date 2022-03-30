FROM python

RUN apt-get update 
RUN apt-get install -y postgresql postgresql-contrib

USER postgres
RUN /etc/init.d/postgresql start && \
    createdb pandas && \
    psql -c "ALTER USER postgres PASSWORD 'postgres';"

USER root
RUN python -m pip install --upgrade pip
RUN python -m pip install \
    cython \
    hypothesis \
    numpy \
    psycopg2 \
    pytest \
    pytest-asyncio \
    python-dateutil \
    pytz \
    sqlalchemy

WORKDIR /pandas


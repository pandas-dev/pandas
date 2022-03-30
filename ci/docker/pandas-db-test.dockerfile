FROM python

RUN apt-get update 
RUN apt-get install -y postgresql postgresql-contrib sqlite3

USER postgres
RUN /etc/init.d/postgresql start && \
    createdb pandas

RUN python -m pip install --upgrade pip
RUN python -m pip install \
    cython \
    numpy \
    pytest \
    python-dateutil \
    pytz
RUN python -m pip install psycopg2

WORKDIR /pandas

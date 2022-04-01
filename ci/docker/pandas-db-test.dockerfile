FROM python

RUN apt-get update

RUN python -m pip install --upgrade pip
RUN python -m pip install \
    cython \
    hypothesis \
    numpy \
    pymysql \
    psycopg2 \
    pytest \
    pytest-asyncio \
    python-dateutil \
    pytz \
    sqlalchemy

WORKDIR /pandas


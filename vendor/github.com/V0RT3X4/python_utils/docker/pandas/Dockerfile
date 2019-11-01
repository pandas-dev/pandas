FROM python:3.7-alpine
ARG PANDAS_VERSION=0.24.1

RUN   apk add libstdc++ \
  &&  apk add --no-cache --virtual .build-deps \
        g++ \
        gcc \
        gfortran \
        build-base \
        wget \
        freetype-dev \
        libpng-dev \
        openblas-dev \
        postgresql-dev \
        musl-dev \
  &&  ln -s /usr/include/locale.h /usr/include/xlocale.h \
  &&  pip install wheel \
  &&  pip install --no-cache \
        numpy==1.15.1 \
        scipy \
        psycopg2-binary \
        sqlalchemy \
  &&  pip install --no-cache \
        pandas==${PANDAS_VERSION} \
  &&  apk del .build-deps

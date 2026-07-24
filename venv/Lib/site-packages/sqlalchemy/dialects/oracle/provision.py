# dialects/oracle/provision.py
# Copyright (C) 2005-2026 the SQLAlchemy authors and contributors
# <see AUTHORS file>
#
# This module is part of SQLAlchemy and is released under
# the MIT License: https://www.opensource.org/licenses/mit-license.php
# mypy: ignore-errors

import time

from ... import create_engine
from ... import exc
from ... import inspect
from ...engine import url as sa_url
from ...testing.provision import configure_follower
from ...testing.provision import create_db
from ...testing.provision import drop_all_schema_objects_post_tables
from ...testing.provision import drop_all_schema_objects_pre_tables
from ...testing.provision import drop_db
from ...testing.provision import follower_url_from_main
from ...testing.provision import generate_driver_url
from ...testing.provision import is_preferred_driver
from ...testing.provision import log
from ...testing.provision import post_configure_engine
from ...testing.provision import post_configure_testing_engine
from ...testing.provision import run_reap_dbs
from ...testing.provision import set_default_schema_on_connection
from ...testing.provision import stop_test_class_outside_fixtures
from ...testing.provision import temp_table_keyword_args
from ...testing.provision import update_db_opts
from ...testing.warnings import warn_test_suite


@generate_driver_url.for_db("oracle")
def _oracle_generate_driver_url(url, driver, query_str):

    backend = url.get_backend_name()

    new_url = url.set(
        drivername="%s+%s" % (backend, driver),
    )

    # use oracledb's retry feature, which is essential for oracle 23c
    # which otherwise frequently rejects connections under load
    # for cx_oracle we have a connect event instead
    if driver in ("oracledb", "oracledb_async"):
        # oracledb is even nice enough to convert from string to int
        # for these opts, apparently
        new_url = new_url.update_query_pairs(
            [("retry_count", "5"), ("retry_delay", "2")]
        )
    else:
        # remove these params for cx_oracle if we received an
        # already-modified URL
        new_url = new_url.difference_update_query(
            ["retry_count", "retry_delay"]
        )

    try:
        new_url.get_dialect()
    except exc.NoSuchModuleError:
        return None
    else:
        return new_url


@create_db.for_db("oracle")
def _oracle_create_db(cfg, eng, ident):
    # NOTE: make sure you've run "ALTER DATABASE default tablespace users" or
    # similar, so that the default tablespace is not "system"; reflection will
    # fail otherwise
    with eng.begin() as conn:
        conn.exec_driver_sql("create user %s identified by xe" % ident)
        conn.exec_driver_sql("create user %s_ts1 identified by xe" % ident)
        conn.exec_driver_sql("create user %s_ts2 identified by xe" % ident)
        conn.exec_driver_sql("grant dba to %s" % (ident,))
        conn.exec_driver_sql("grant unlimited tablespace to %s" % ident)
        conn.exec_driver_sql("grant unlimited tablespace to %s_ts1" % ident)
        conn.exec_driver_sql("grant unlimited tablespace to %s_ts2" % ident)
        # these are needed to create materialized views
        conn.exec_driver_sql("grant create table to %s" % ident)
        conn.exec_driver_sql("grant create table to %s_ts1" % ident)
        conn.exec_driver_sql("grant create table to %s_ts2" % ident)


@configure_follower.for_db("oracle")
def _oracle_configure_follower(config, ident):
    config.test_schema = "%s_ts1" % ident
    config.test_schema_2 = "%s_ts2" % ident


def _ora_drop_ignore(conn, dbname):
    try:
        conn.exec_driver_sql("drop user %s cascade" % dbname)
        log.info("Reaped db: %s", dbname)
        return True
    except exc.DatabaseError as err:
        log.warning("couldn't drop db: %s", err)
        return False


@drop_all_schema_objects_pre_tables.for_db("oracle")
def _ora_drop_all_schema_objects_pre_tables(cfg, eng):
    _purge_recyclebin(eng)
    _purge_recyclebin(eng, cfg.test_schema)


@drop_all_schema_objects_post_tables.for_db("oracle")
def _ora_drop_all_schema_objects_post_tables(cfg, eng):
    with eng.begin() as conn:
        for syn in conn.dialect._get_synonyms(conn, None, None, None):
            conn.exec_driver_sql(f"drop synonym {syn['synonym_name']}")

        for syn in conn.dialect._get_synonyms(
            conn, cfg.test_schema, None, None
        ):
            conn.exec_driver_sql(
                f"drop synonym {cfg.test_schema}.{syn['synonym_name']}"
            )

        for tmp_table in inspect(conn).get_temp_table_names():
            conn.exec_driver_sql(f"drop table {tmp_table}")


@drop_db.for_db("oracle")
def _oracle_drop_db(cfg, eng, ident):
    with eng.begin() as conn:
        # cx_Oracle seems to occasionally leak open connections when a large
        # suite it run, even if we confirm we have zero references to
        # connection objects.
        # while there is a "kill session" command in Oracle Database,
        # it unfortunately does not release the connection sufficiently.
        _ora_drop_ignore(conn, ident)
        _ora_drop_ignore(conn, "%s_ts1" % ident)
        _ora_drop_ignore(conn, "%s_ts2" % ident)


@stop_test_class_outside_fixtures.for_db("oracle")
def _ora_stop_test_class_outside_fixtures(config, db, cls):
    try:
        _purge_recyclebin(db)
    except exc.DatabaseError as err:
        log.warning("purge recyclebin command failed: %s", err)


def _purge_recyclebin(eng, schema=None):
    with eng.begin() as conn:
        if schema is None:
            # run magic command to get rid of identity sequences
            # https://floo.bar/2019/11/29/drop-the-underlying-sequence-of-an-identity-column/  # noqa: E501
            conn.exec_driver_sql("purge recyclebin")
        else:
            # per user: https://community.oracle.com/tech/developers/discussion/2255402/how-to-clear-dba-recyclebin-for-a-particular-user  # noqa: E501
            for owner, object_name, type_ in conn.exec_driver_sql(
                "select owner, object_name,type from "
                "dba_recyclebin where owner=:schema and type='TABLE'",
                {"schema": conn.dialect.denormalize_name(schema)},
            ).all():
                conn.exec_driver_sql(f'purge {type_} {owner}."{object_name}"')


@is_preferred_driver.for_db("oracle")
def _oracle_is_preferred_driver(cfg, engine):
    """establish oracledb as the preferred driver to use for tests, even
    though cx_Oracle is still the "default" driver"""

    return engine.dialect.driver == "oracledb" and not engine.dialect.is_async


def _connect_with_retry(dialect, conn_rec, cargs, cparams):
    assert dialect.driver == "cx_oracle"

    def _is_couldnt_connect(err):
        return "DPY-6005" in str(err) or "ORA-12516" in str(err)

    err_ = None
    for _ in range(5):
        try:
            return dialect.loaded_dbapi.connect(*cargs, **cparams)
        except (
            dialect.loaded_dbapi.DatabaseError,
            dialect.loaded_dbapi.OperationalError,
        ) as err:
            err_ = err
            if _is_couldnt_connect(err):
                warn_test_suite("Oracle database reconnecting...")
                time.sleep(2)
                continue
            else:
                raise
    if err_ is not None:
        raise Exception("connect failed after five attempts") from err_


@post_configure_testing_engine.for_db("oracle")
def _oracle_post_configure_testing_engine(url, engine, options, scope):
    from ... import event

    if engine.dialect.driver == "cx_oracle":
        event.listen(engine, "do_connect", _connect_with_retry)


@post_configure_engine.for_db("oracle")
def _oracle_post_configure_engine(url, engine, follower_ident):

    from ... import event

    @event.listens_for(engine, "checkin")
    def checkin(dbapi_connection, connection_record):
        # this was meant to work around this issue:
        # https://github.com/oracle/python-cx_Oracle/issues/530
        # invalidate oracle connections that had 2pc set up
        # however things are too complex with some of the 2pc tests,
        # so just block cx_oracle from being used in 2pc tests (use oracledb
        # instead)
        # if "cx_oracle_xid" in connection_record.info:
        #    connection_record.invalidate()

        # clear statement cache on all connections that were used
        # https://github.com/oracle/python-cx_Oracle/issues/519
        # TODO: oracledb claims to have this feature built in somehow,
        # see if that's in use and/or if it needs to be enabled
        # (or if this doesn't even apply to the newer oracle's we're using)
        try:
            sc = dbapi_connection.stmtcachesize
        except:
            # connection closed
            pass
        else:
            dbapi_connection.stmtcachesize = 0
            dbapi_connection.stmtcachesize = sc


@run_reap_dbs.for_db("oracle")
def _reap_oracle_dbs(url, idents):
    log.info("db reaper connecting to %r", url)
    eng = create_engine(url)
    with eng.begin() as conn:
        log.info("identifiers in file: %s", ", ".join(idents))

        to_reap = conn.exec_driver_sql(
            "select u.username from all_users u where username "
            "like 'TEST_%' and not exists (select username "
            "from v$session where username=u.username)"
        )
        all_names = {username.lower() for (username,) in to_reap}
        to_drop = set()
        for name in all_names:
            if name.endswith("_ts1") or name.endswith("_ts2"):
                continue
            elif name in idents:
                to_drop.add(name)
                if "%s_ts1" % name in all_names:
                    to_drop.add("%s_ts1" % name)
                if "%s_ts2" % name in all_names:
                    to_drop.add("%s_ts2" % name)

        dropped = total = 0
        for total, username in enumerate(to_drop, 1):
            if _ora_drop_ignore(conn, username):
                dropped += 1
        log.info(
            "Dropped %d out of %d stale databases detected", dropped, total
        )


@follower_url_from_main.for_db("oracle")
def _oracle_follower_url_from_main(url, ident):
    url = sa_url.make_url(url)
    return url.set(username=ident, password="xe")


@temp_table_keyword_args.for_db("oracle")
def _oracle_temp_table_keyword_args(cfg, eng):
    return {
        "prefixes": ["GLOBAL TEMPORARY"],
        "oracle_on_commit": "PRESERVE ROWS",
    }


@set_default_schema_on_connection.for_db("oracle")
def _oracle_set_default_schema_on_connection(
    cfg, dbapi_connection, schema_name
):
    cursor = dbapi_connection.cursor()
    cursor.execute("ALTER SESSION SET CURRENT_SCHEMA=%s" % schema_name)
    cursor.close()


@update_db_opts.for_db("oracle")
def _update_db_opts(db_url, db_opts, options):
    """Set database options (db_opts) for a test database that we created."""
    if (
        options.oracledb_thick_mode
        and sa_url.make_url(db_url).get_driver_name() == "oracledb"
    ):
        db_opts["thick_mode"] = True

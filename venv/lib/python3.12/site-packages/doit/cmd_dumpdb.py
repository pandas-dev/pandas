import pprint
import json
import dbm
from dbm import whichdb


from .exceptions import InvalidCommand
from .cmd_base import Command, opt_depfile


def dbm_iter(db):
    # try dictionary interface - used in dumbdb
    try:
        return db.items()
    except AttributeError:  # pragma: no cover
        pass

    # try firstkey/nextkey - ok for py3 dbm.gnu
    try:  # pragma: no cover
        db.firstkey
        def iter_gdbm(db):
            k = db.firstkey()
            while k is not None:
                yield k, db[k]
                k = db.nextkey(k)
        return iter_gdbm(db)
    except Exception:  # pragma: no cover
        raise InvalidCommand("It seems your DB backend doesn't support "
                             "iterating through all elements")


class DumpDB(Command):
    """dump dependency DB"""
    doc_purpose = 'dump dependency DB'
    doc_usage = ''
    doc_description = None

    cmd_options = (opt_depfile,)

    def execute(self, opt_values, pos_args):
        dep_file = opt_values['dep_file']
        db_type = whichdb(dep_file)
        print("DBM type is '%s'" % db_type)
        if db_type in ('dbm', 'dbm.ndbm'):  # pragma: no cover
            raise InvalidCommand('ndbm does not support iteration of elements')
        data = dbm.open(dep_file)
        for key, value_str in dbm_iter(data):
            value_dict = json.loads(value_str.decode('utf-8'))
            value_fmt = pprint.pformat(value_dict, indent=4, width=100)
            print("{key} -> {value}".format(key=key, value=value_fmt))

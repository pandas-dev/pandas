# Python implementation of the MySQL client-server protocol
# http://dev.mysql.com/doc/internals/en/client-server-protocol.html
# Error codes:
# https://dev.mysql.com/doc/refman/5.5/en/error-handling.html
import errno
import os
import socket
import struct
import sys
import traceback
import warnings

from . import _auth

from .charset import charset_by_name, charset_by_id
from .constants import CLIENT, COMMAND, CR, ER, FIELD_TYPE, SERVER_STATUS
from . import converters
from .cursors import Cursor
from .optionfile import Parser
from .protocol import (
    dump_packet,
    MysqlPacket,
    FieldDescriptorPacket,
    OKPacketWrapper,
    EOFPacketWrapper,
    LoadLocalPacketWrapper,
)
from . import err, VERSION_STRING

try:
    import ssl

    SSL_ENABLED = True
except ImportError:
    ssl = None
    SSL_ENABLED = False

try:
    import getpass

    DEFAULT_USER = getpass.getuser()
    del getpass
except (ImportError, KeyError):
    # KeyError occurs when there's no entry in OS database for a current user.
    DEFAULT_USER = None

DEBUG = False

TEXT_TYPES = {
    FIELD_TYPE.BIT,
    FIELD_TYPE.BLOB,
    FIELD_TYPE.LONG_BLOB,
    FIELD_TYPE.MEDIUM_BLOB,
    FIELD_TYPE.STRING,
    FIELD_TYPE.TINY_BLOB,
    FIELD_TYPE.VAR_STRING,
    FIELD_TYPE.VARCHAR,
    FIELD_TYPE.GEOMETRY,
}


DEFAULT_CHARSET = "utf8mb4"

MAX_PACKET_LEN = 2**24 - 1


def _pack_int24(n):
    return struct.pack("<I", n)[:3]


# https://dev.mysql.com/doc/internals/en/integer.html#packet-Protocol::LengthEncodedInteger
def _lenenc_int(i):
    if i < 0:
        raise ValueError(
            "Encoding %d is less than 0 - no representation in LengthEncodedInteger" % i
        )
    elif i < 0xFB:
        return bytes([i])
    elif i < (1 << 16):
        return b"\xfc" + struct.pack("<H", i)
    elif i < (1 << 24):
        return b"\xfd" + struct.pack("<I", i)[:3]
    elif i < (1 << 64):
        return b"\xfe" + struct.pack("<Q", i)
    else:
        raise ValueError(
            f"Encoding {i:x} is larger than {1 << 64:x} - no representation in LengthEncodedInteger"
        )


class Connection:
    """
    Representation of a socket with a mysql server.

    The proper way to get an instance of this class is to call
    connect().

    Establish a connection to the MySQL database. Accepts several
    arguments:

    :param host: Host where the database server is located.
    :param user: Username to log in as.
    :param password: Password to use.
    :param database: Database to use, None to not use a particular one.
    :param port: MySQL port to use, default is usually OK. (default: 3306)
    :param bind_address: When the client has multiple network interfaces, specify
        the interface from which to connect to the host. Argument can be
        a hostname or an IP address.
    :param unix_socket: Use a unix socket rather than TCP/IP.
    :param read_timeout: The timeout for reading from the connection in seconds.
        (default: None - no timeout)
    :param write_timeout: The timeout for writing to the connection in seconds.
        (default: None - no timeout)
    :param str charset: Charset to use.
    :param str collation: Collation name to use.
    :param sql_mode: Default SQL_MODE to use.
    :param read_default_file:
        Specifies  my.cnf file to read these parameters from under the [client] section.
    :param conv:
        Conversion dictionary to use instead of the default one.
        This is used to provide custom marshalling and unmarshalling of types.
        See converters.
    :param use_unicode:
        Whether or not to default to unicode strings.
        This option defaults to true.
    :param client_flag: Custom flags to send to MySQL. Find potential values in constants.CLIENT.
    :param cursorclass: Custom cursor class to use.
    :param init_command: Initial SQL statement to run when connection is established.
    :param connect_timeout: The timeout for connecting to the database in seconds.
        (default: 10, min: 1, max: 31536000)
    :param ssl: A dict of arguments similar to mysql_ssl_set()'s parameters or an ssl.SSLContext.
    :param ssl_ca: Path to the file that contains a PEM-formatted CA certificate.
    :param ssl_cert: Path to the file that contains a PEM-formatted client certificate.
    :param ssl_disabled: A boolean value that disables usage of TLS.
    :param ssl_key: Path to the file that contains a PEM-formatted private key for
        the client certificate.
    :param ssl_key_password: The password for the client certificate private key.
    :param ssl_verify_cert: Set to true to check the server certificate's validity.
    :param ssl_verify_identity: Set to true to check the server's identity.
    :param read_default_group: Group to read from in the configuration file.
    :param autocommit: Autocommit mode. None means use server default. (default: False)
    :param local_infile: Boolean to enable the use of LOAD DATA LOCAL command. (default: False)
    :param max_allowed_packet: Max size of packet sent to server in bytes. (default: 16MB)
        Only used to limit size of "LOAD LOCAL INFILE" data packet smaller than default (16KB).
    :param defer_connect: Don't explicitly connect on construction - wait for connect call.
        (default: False)
    :param auth_plugin_map: A dict of plugin names to a class that processes that plugin.
        The class will take the Connection object as the argument to the constructor.
        The class needs an authenticate method taking an authentication packet as
        an argument.  For the dialog plugin, a prompt(echo, prompt) method can be used
        (if no authenticate method) for returning a string from the user. (experimental)
    :param server_public_key: SHA256 authentication plugin public key value. (default: None)
    :param binary_prefix: Add _binary prefix on bytes and bytearray. (default: False)
    :param compress: Not supported.
    :param named_pipe: Not supported.
    :param db: **DEPRECATED** Alias for database.
    :param passwd: **DEPRECATED** Alias for password.

    See `Connection <https://www.python.org/dev/peps/pep-0249/#connection-objects>`_ in the
    specification.
    """

    _sock = None
    _auth_plugin_name = ""
    _closed = False
    _secure = False

    def __init__(
        self,
        *,
        user=None,  # The first four arguments is based on DB-API 2.0 recommendation.
        password="",
        host=None,
        database=None,
        unix_socket=None,
        port=0,
        charset="",
        collation=None,
        sql_mode=None,
        read_default_file=None,
        conv=None,
        use_unicode=True,
        client_flag=0,
        cursorclass=Cursor,
        init_command=None,
        connect_timeout=10,
        read_default_group=None,
        autocommit=False,
        local_infile=False,
        max_allowed_packet=16 * 1024 * 1024,
        defer_connect=False,
        auth_plugin_map=None,
        read_timeout=None,
        write_timeout=None,
        bind_address=None,
        binary_prefix=False,
        program_name=None,
        server_public_key=None,
        ssl=None,
        ssl_ca=None,
        ssl_cert=None,
        ssl_disabled=None,
        ssl_key=None,
        ssl_key_password=None,
        ssl_verify_cert=None,
        ssl_verify_identity=None,
        compress=None,  # not supported
        named_pipe=None,  # not supported
        passwd=None,  # deprecated
        db=None,  # deprecated
    ):
        if db is not None and database is None:
            # We will raise warning in 2022 or later.
            # See https://github.com/PyMySQL/PyMySQL/issues/939
            # warnings.warn("'db' is deprecated, use 'database'", DeprecationWarning, 3)
            database = db
        if passwd is not None and not password:
            # We will raise warning in 2022 or later.
            # See https://github.com/PyMySQL/PyMySQL/issues/939
            # warnings.warn(
            #    "'passwd' is deprecated, use 'password'", DeprecationWarning, 3
            # )
            password = passwd

        if compress or named_pipe:
            raise NotImplementedError(
                "compress and named_pipe arguments are not supported"
            )

        self._local_infile = bool(local_infile)
        if self._local_infile:
            client_flag |= CLIENT.LOCAL_FILES

        if read_default_group and not read_default_file:
            if sys.platform.startswith("win"):
                read_default_file = "c:\\my.ini"
            else:
                read_default_file = "/etc/my.cnf"

        if read_default_file:
            if not read_default_group:
                read_default_group = "client"

            cfg = Parser()
            cfg.read(os.path.expanduser(read_default_file))

            def _config(key, arg):
                if arg:
                    return arg
                try:
                    return cfg.get(read_default_group, key)
                except Exception:
                    return arg

            user = _config("user", user)
            password = _config("password", password)
            host = _config("host", host)
            database = _config("database", database)
            unix_socket = _config("socket", unix_socket)
            port = int(_config("port", port))
            bind_address = _config("bind-address", bind_address)
            charset = _config("default-character-set", charset)
            if not ssl:
                ssl = {}
            if isinstance(ssl, dict):
                for key in ["ca", "capath", "cert", "key", "password", "cipher"]:
                    value = _config("ssl-" + key, ssl.get(key))
                    if value:
                        ssl[key] = value

        self.ssl = False
        if not ssl_disabled:
            if ssl_ca or ssl_cert or ssl_key or ssl_verify_cert or ssl_verify_identity:
                ssl = {
                    "ca": ssl_ca,
                    "check_hostname": bool(ssl_verify_identity),
                    "verify_mode": ssl_verify_cert
                    if ssl_verify_cert is not None
                    else False,
                }
                if ssl_cert is not None:
                    ssl["cert"] = ssl_cert
                if ssl_key is not None:
                    ssl["key"] = ssl_key
                if ssl_key_password is not None:
                    ssl["password"] = ssl_key_password
            if ssl:
                if not SSL_ENABLED:
                    raise NotImplementedError("ssl module not found")
                self.ssl = True
                client_flag |= CLIENT.SSL
                self.ctx = self._create_ssl_ctx(ssl)

        self.host = host or "localhost"
        self.port = port or 3306
        if type(self.port) is not int:
            raise ValueError("port should be of type int")
        self.user = user or DEFAULT_USER
        self.password = password or b""
        if isinstance(self.password, str):
            self.password = self.password.encode("latin1")
        self.db = database
        self.unix_socket = unix_socket
        self.bind_address = bind_address
        if not (0 < connect_timeout <= 31536000):
            raise ValueError("connect_timeout should be >0 and <=31536000")
        self.connect_timeout = connect_timeout or None
        if read_timeout is not None and read_timeout <= 0:
            raise ValueError("read_timeout should be > 0")
        self._read_timeout = read_timeout
        if write_timeout is not None and write_timeout <= 0:
            raise ValueError("write_timeout should be > 0")
        self._write_timeout = write_timeout

        self.charset = charset or DEFAULT_CHARSET
        self.collation = collation
        self.use_unicode = use_unicode

        self.encoding = charset_by_name(self.charset).encoding

        client_flag |= CLIENT.CAPABILITIES
        if self.db:
            client_flag |= CLIENT.CONNECT_WITH_DB

        self.client_flag = client_flag

        self.cursorclass = cursorclass

        self._result = None
        self._affected_rows = 0
        self.host_info = "Not connected"

        # specified autocommit mode. None means use server default.
        self.autocommit_mode = autocommit

        if conv is None:
            conv = converters.conversions

        # Need for MySQLdb compatibility.
        self.encoders = {k: v for (k, v) in conv.items() if type(k) is not int}
        self.decoders = {k: v for (k, v) in conv.items() if type(k) is int}
        self.sql_mode = sql_mode
        self.init_command = init_command
        self.max_allowed_packet = max_allowed_packet
        self._auth_plugin_map = auth_plugin_map or {}
        self._binary_prefix = binary_prefix
        self.server_public_key = server_public_key

        self._connect_attrs = {
            "_client_name": "pymysql",
            "_client_version": VERSION_STRING,
            "_pid": str(os.getpid()),
        }

        if program_name:
            self._connect_attrs["program_name"] = program_name

        if defer_connect:
            self._sock = None
        else:
            self.connect()

    def __enter__(self):
        return self

    def __exit__(self, *exc_info):
        del exc_info
        self.close()

    def _create_ssl_ctx(self, sslp):
        if isinstance(sslp, ssl.SSLContext):
            return sslp
        ca = sslp.get("ca")
        capath = sslp.get("capath")
        hasnoca = ca is None and capath is None
        ctx = ssl.create_default_context(cafile=ca, capath=capath)
        ctx.check_hostname = not hasnoca and sslp.get("check_hostname", True)
        verify_mode_value = sslp.get("verify_mode")
        if verify_mode_value is None:
            ctx.verify_mode = ssl.CERT_NONE if hasnoca else ssl.CERT_REQUIRED
        elif isinstance(verify_mode_value, bool):
            ctx.verify_mode = ssl.CERT_REQUIRED if verify_mode_value else ssl.CERT_NONE
        else:
            if isinstance(verify_mode_value, str):
                verify_mode_value = verify_mode_value.lower()
            if verify_mode_value in ("none", "0", "false", "no"):
                ctx.verify_mode = ssl.CERT_NONE
            elif verify_mode_value == "optional":
                ctx.verify_mode = ssl.CERT_OPTIONAL
            elif verify_mode_value in ("required", "1", "true", "yes"):
                ctx.verify_mode = ssl.CERT_REQUIRED
            else:
                ctx.verify_mode = ssl.CERT_NONE if hasnoca else ssl.CERT_REQUIRED
        if "cert" in sslp:
            ctx.load_cert_chain(
                sslp["cert"], keyfile=sslp.get("key"), password=sslp.get("password")
            )
        if "cipher" in sslp:
            ctx.set_ciphers(sslp["cipher"])
        ctx.options |= ssl.OP_NO_SSLv2
        ctx.options |= ssl.OP_NO_SSLv3
        return ctx

    def close(self):
        """
        Send the quit message and close the socket.

        See `Connection.close() <https://www.python.org/dev/peps/pep-0249/#Connection.close>`_
        in the specification.

        :raise Error: If the connection is already closed.
        """
        if self._closed:
            raise err.Error("Already closed")
        self._closed = True
        if self._sock is None:
            return
        send_data = struct.pack("<iB", 1, COMMAND.COM_QUIT)
        try:
            self._write_bytes(send_data)
        except Exception:
            pass
        finally:
            self._force_close()

    @property
    def open(self):
        """Return True if the connection is open."""
        return self._sock is not None

    def _force_close(self):
        """Close connection without QUIT message."""
        if self._sock:
            try:
                self._sock.close()
            except:  # noqa
                pass
        self._sock = None
        self._rfile = None

    __del__ = _force_close

    def autocommit(self, value):
        self.autocommit_mode = bool(value)
        current = self.get_autocommit()
        if value != current:
            self._send_autocommit_mode()

    def get_autocommit(self):
        return bool(self.server_status & SERVER_STATUS.SERVER_STATUS_AUTOCOMMIT)

    def _read_ok_packet(self):
        pkt = self._read_packet()
        if not pkt.is_ok_packet():
            raise err.OperationalError(
                CR.CR_COMMANDS_OUT_OF_SYNC,
                "Command Out of Sync",
            )
        ok = OKPacketWrapper(pkt)
        self.server_status = ok.server_status
        return ok

    def _send_autocommit_mode(self):
        """Set whether or not to commit after every execute()."""
        self._execute_command(
            COMMAND.COM_QUERY, "SET AUTOCOMMIT = %s" % self.escape(self.autocommit_mode)
        )
        self._read_ok_packet()

    def begin(self):
        """Begin transaction."""
        self._execute_command(COMMAND.COM_QUERY, "BEGIN")
        self._read_ok_packet()

    def commit(self):
        """
        Commit changes to stable storage.

        See `Connection.commit() <https://www.python.org/dev/peps/pep-0249/#commit>`_
        in the specification.
        """
        self._execute_command(COMMAND.COM_QUERY, "COMMIT")
        self._read_ok_packet()

    def rollback(self):
        """
        Roll back the current transaction.

        See `Connection.rollback() <https://www.python.org/dev/peps/pep-0249/#rollback>`_
        in the specification.
        """
        self._execute_command(COMMAND.COM_QUERY, "ROLLBACK")
        self._read_ok_packet()

    def show_warnings(self):
        """Send the "SHOW WARNINGS" SQL command."""
        self._execute_command(COMMAND.COM_QUERY, "SHOW WARNINGS")
        result = MySQLResult(self)
        result.read()
        return result.rows

    def select_db(self, db):
        """
        Set current db.

        :param db: The name of the db.
        """
        self._execute_command(COMMAND.COM_INIT_DB, db)
        self._read_ok_packet()

    def escape(self, obj, mapping=None):
        """Escape whatever value is passed.

        Non-standard, for internal use; do not use this in your applications.
        """
        if isinstance(obj, str):
            return "'" + self.escape_string(obj) + "'"
        if isinstance(obj, (bytes, bytearray)):
            ret = self._quote_bytes(obj)
            if self._binary_prefix:
                ret = "_binary" + ret
            return ret
        return converters.escape_item(obj, self.charset, mapping=mapping)

    def literal(self, obj):
        """Alias for escape().

        Non-standard, for internal use; do not use this in your applications.
        """
        return self.escape(obj, self.encoders)

    def escape_string(self, s):
        if self.server_status & SERVER_STATUS.SERVER_STATUS_NO_BACKSLASH_ESCAPES:
            return s.replace("'", "''")
        return converters.escape_string(s)

    def _quote_bytes(self, s):
        if self.server_status & SERVER_STATUS.SERVER_STATUS_NO_BACKSLASH_ESCAPES:
            return "'{}'".format(
                s.replace(b"'", b"''").decode("ascii", "surrogateescape")
            )
        return converters.escape_bytes(s)

    def cursor(self, cursor=None):
        """
        Create a new cursor to execute queries with.

        :param cursor: The type of cursor to create. None means use Cursor.
        :type cursor: :py:class:`Cursor`, :py:class:`SSCursor`, :py:class:`DictCursor`,
            or :py:class:`SSDictCursor`.
        """
        if cursor:
            return cursor(self)
        return self.cursorclass(self)

    # The following methods are INTERNAL USE ONLY (called from Cursor)
    def query(self, sql, unbuffered=False):
        # if DEBUG:
        #     print("DEBUG: sending query:", sql)
        if isinstance(sql, str):
            sql = sql.encode(self.encoding, "surrogateescape")
        self._execute_command(COMMAND.COM_QUERY, sql)
        self._affected_rows = self._read_query_result(unbuffered=unbuffered)
        return self._affected_rows

    def next_result(self, unbuffered=False):
        self._affected_rows = self._read_query_result(unbuffered=unbuffered)
        return self._affected_rows

    def affected_rows(self):
        return self._affected_rows

    def kill(self, thread_id):
        arg = struct.pack("<I", thread_id)
        self._execute_command(COMMAND.COM_PROCESS_KILL, arg)
        return self._read_ok_packet()

    def ping(self, reconnect=True):
        """
        Check if the server is alive.

        :param reconnect: If the connection is closed, reconnect.
        :type reconnect: boolean

        :raise Error: If the connection is closed and reconnect=False.
        """
        if self._sock is None:
            if reconnect:
                self.connect()
                reconnect = False
            else:
                raise err.Error("Already closed")
        try:
            self._execute_command(COMMAND.COM_PING, "")
            self._read_ok_packet()
        except Exception:
            if reconnect:
                self.connect()
                self.ping(False)
            else:
                raise

    def set_charset(self, charset):
        """Deprecated. Use set_character_set() instead."""
        # This function has been implemented in old PyMySQL.
        # But this name is different from MySQLdb.
        # So we keep this function for compatibility and add
        # new set_character_set() function.
        self.set_character_set(charset)

    def set_character_set(self, charset, collation=None):
        """
        Set charaset (and collation)

        Send "SET NAMES charset [COLLATE collation]" query.
        Update Connection.encoding based on charset.
        """
        # Make sure charset is supported.
        encoding = charset_by_name(charset).encoding

        if collation:
            query = f"SET NAMES {charset} COLLATE {collation}"
        else:
            query = f"SET NAMES {charset}"
        self._execute_command(COMMAND.COM_QUERY, query)
        self._read_packet()
        self.charset = charset
        self.encoding = encoding
        self.collation = collation

    def connect(self, sock=None):
        self._closed = False
        try:
            if sock is None:
                if self.unix_socket:
                    sock = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
                    sock.settimeout(self.connect_timeout)
                    sock.connect(self.unix_socket)
                    self.host_info = "Localhost via UNIX socket"
                    self._secure = True
                    if DEBUG:
                        print("connected using unix_socket")
                else:
                    kwargs = {}
                    if self.bind_address is not None:
                        kwargs["source_address"] = (self.bind_address, 0)
                    while True:
                        try:
                            sock = socket.create_connection(
                                (self.host, self.port), self.connect_timeout, **kwargs
                            )
                            break
                        except OSError as e:
                            if e.errno == errno.EINTR:
                                continue
                            raise
                    self.host_info = "socket %s:%d" % (self.host, self.port)
                    if DEBUG:
                        print("connected using socket")
                    sock.setsockopt(socket.IPPROTO_TCP, socket.TCP_NODELAY, 1)
                    sock.setsockopt(socket.SOL_SOCKET, socket.SO_KEEPALIVE, 1)
                sock.settimeout(None)

            self._sock = sock
            self._rfile = sock.makefile("rb")
            self._next_seq_id = 0

            self._get_server_information()
            self._request_authentication()

            # Send "SET NAMES" query on init for:
            # - Ensure charaset (and collation) is set to the server.
            #   - collation_id in handshake packet may be ignored.
            # - If collation is not specified, we don't know what is server's
            #   default collation for the charset. For example, default collation
            #   of utf8mb4 is:
            #   - MySQL 5.7, MariaDB 10.x: utf8mb4_general_ci
            #   - MySQL 8.0: utf8mb4_0900_ai_ci
            #
            # Reference:
            # - https://github.com/PyMySQL/PyMySQL/issues/1092
            # - https://github.com/wagtail/wagtail/issues/9477
            # - https://zenn.dev/methane/articles/2023-mysql-collation (Japanese)
            self.set_character_set(self.charset, self.collation)

            if self.sql_mode is not None:
                c = self.cursor()
                c.execute("SET sql_mode=%s", (self.sql_mode,))
                c.close()

            if self.init_command is not None:
                c = self.cursor()
                c.execute(self.init_command)
                c.close()

            if self.autocommit_mode is not None:
                self.autocommit(self.autocommit_mode)
        except BaseException as e:
            self._rfile = None
            if sock is not None:
                try:
                    sock.close()
                except:  # noqa
                    pass

            if isinstance(e, (OSError, IOError)):
                exc = err.OperationalError(
                    CR.CR_CONN_HOST_ERROR,
                    f"Can't connect to MySQL server on {self.host!r} ({e})",
                )
                # Keep original exception and traceback to investigate error.
                exc.original_exception = e
                exc.traceback = traceback.format_exc()
                if DEBUG:
                    print(exc.traceback)
                raise exc

            # If e is neither DatabaseError or IOError, It's a bug.
            # But raising AssertionError hides original error.
            # So just reraise it.
            raise

    def write_packet(self, payload):
        """Writes an entire "mysql packet" in its entirety to the network
        adding its length and sequence number.
        """
        # Internal note: when you build packet manually and calls _write_bytes()
        # directly, you should set self._next_seq_id properly.
        data = _pack_int24(len(payload)) + bytes([self._next_seq_id]) + payload
        if DEBUG:
            dump_packet(data)
        self._write_bytes(data)
        self._next_seq_id = (self._next_seq_id + 1) % 256

    def _read_packet(self, packet_type=MysqlPacket):
        """Read an entire "mysql packet" in its entirety from the network
        and return a MysqlPacket type that represents the results.

        :raise OperationalError: If the connection to the MySQL server is lost.
        :raise InternalError: If the packet sequence number is wrong.
        """
        buff = bytearray()
        while True:
            packet_header = self._read_bytes(4)
            # if DEBUG: dump_packet(packet_header)

            btrl, btrh, packet_number = struct.unpack("<HBB", packet_header)
            bytes_to_read = btrl + (btrh << 16)
            if packet_number != self._next_seq_id:
                self._force_close()
                if packet_number == 0:
                    # MariaDB sends error packet with seqno==0 when shutdown
                    raise err.OperationalError(
                        CR.CR_SERVER_LOST,
                        "Lost connection to MySQL server during query",
                    )
                raise err.InternalError(
                    "Packet sequence number wrong - got %d expected %d"
                    % (packet_number, self._next_seq_id)
                )
            self._next_seq_id = (self._next_seq_id + 1) % 256

            recv_data = self._read_bytes(bytes_to_read)
            if DEBUG:
                dump_packet(recv_data)
            buff += recv_data
            # https://dev.mysql.com/doc/internals/en/sending-more-than-16mbyte.html
            if bytes_to_read < MAX_PACKET_LEN:
                break

        packet = packet_type(bytes(buff), self.encoding)
        if packet.is_error_packet():
            if self._result is not None and self._result.unbuffered_active is True:
                self._result.unbuffered_active = False
            packet.raise_for_error()
        return packet

    def _read_bytes(self, num_bytes):
        self._sock.settimeout(self._read_timeout)
        while True:
            try:
                data = self._rfile.read(num_bytes)
                break
            except OSError as e:
                if e.errno == errno.EINTR:
                    continue
                self._force_close()
                raise err.OperationalError(
                    CR.CR_SERVER_LOST,
                    f"Lost connection to MySQL server during query ({e})",
                )
            except BaseException:
                # Don't convert unknown exception to MySQLError.
                self._force_close()
                raise
        if len(data) < num_bytes:
            self._force_close()
            raise err.OperationalError(
                CR.CR_SERVER_LOST, "Lost connection to MySQL server during query"
            )
        return data

    def _write_bytes(self, data):
        self._sock.settimeout(self._write_timeout)
        try:
            self._sock.sendall(data)
        except OSError as e:
            self._force_close()
            raise err.OperationalError(
                CR.CR_SERVER_GONE_ERROR, f"MySQL server has gone away ({e!r})"
            )

    def _read_query_result(self, unbuffered=False):
        self._result = None
        if unbuffered:
            try:
                result = MySQLResult(self)
                result.init_unbuffered_query()
            except:
                result.unbuffered_active = False
                result.connection = None
                raise
        else:
            result = MySQLResult(self)
            result.read()
        self._result = result
        if result.server_status is not None:
            self.server_status = result.server_status
        return result.affected_rows

    def insert_id(self):
        if self._result:
            return self._result.insert_id
        else:
            return 0

    def _execute_command(self, command, sql):
        """
        :raise InterfaceError: If the connection is closed.
        :raise ValueError: If no username was specified.
        """
        if not self._sock:
            raise err.InterfaceError(0, "")

        # If the last query was unbuffered, make sure it finishes before
        # sending new commands
        if self._result is not None:
            if self._result.unbuffered_active:
                warnings.warn("Previous unbuffered result was left incomplete")
                self._result._finish_unbuffered_query()
            while self._result.has_next:
                self.next_result()
            self._result = None

        if isinstance(sql, str):
            sql = sql.encode(self.encoding)

        packet_size = min(MAX_PACKET_LEN, len(sql) + 1)  # +1 is for command

        # tiny optimization: build first packet manually instead of
        # calling self..write_packet()
        prelude = struct.pack("<iB", packet_size, command)
        packet = prelude + sql[: packet_size - 1]
        self._write_bytes(packet)
        if DEBUG:
            dump_packet(packet)
        self._next_seq_id = 1

        if packet_size < MAX_PACKET_LEN:
            return

        sql = sql[packet_size - 1 :]
        while True:
            packet_size = min(MAX_PACKET_LEN, len(sql))
            self.write_packet(sql[:packet_size])
            sql = sql[packet_size:]
            if not sql and packet_size < MAX_PACKET_LEN:
                break

    def _request_authentication(self):
        # https://dev.mysql.com/doc/internals/en/connection-phase-packets.html#packet-Protocol::HandshakeResponse
        if int(self.server_version.split(".", 1)[0]) >= 5:
            self.client_flag |= CLIENT.MULTI_RESULTS

        if self.user is None:
            raise ValueError("Did not specify a username")

        charset_id = charset_by_name(self.charset).id
        if isinstance(self.user, str):
            self.user = self.user.encode(self.encoding)

        data_init = struct.pack(
            "<iIB23s", self.client_flag, MAX_PACKET_LEN, charset_id, b""
        )

        if self.ssl and self.server_capabilities & CLIENT.SSL:
            self.write_packet(data_init)

            self._sock = self.ctx.wrap_socket(self._sock, server_hostname=self.host)
            self._rfile = self._sock.makefile("rb")
            self._secure = True

        data = data_init + self.user + b"\0"

        authresp = b""
        plugin_name = None

        if self._auth_plugin_name == "":
            plugin_name = b""
            authresp = _auth.scramble_native_password(self.password, self.salt)
        elif self._auth_plugin_name == "mysql_native_password":
            plugin_name = b"mysql_native_password"
            authresp = _auth.scramble_native_password(self.password, self.salt)
        elif self._auth_plugin_name == "caching_sha2_password":
            plugin_name = b"caching_sha2_password"
            if self.password:
                if DEBUG:
                    print("caching_sha2: trying fast path")
                authresp = _auth.scramble_caching_sha2(self.password, self.salt)
            else:
                if DEBUG:
                    print("caching_sha2: empty password")
        elif self._auth_plugin_name == "sha256_password":
            plugin_name = b"sha256_password"
            if self.ssl and self.server_capabilities & CLIENT.SSL:
                authresp = self.password + b"\0"
            elif self.password:
                authresp = b"\1"  # request public key
            else:
                authresp = b"\0"  # empty password

        if self.server_capabilities & CLIENT.PLUGIN_AUTH_LENENC_CLIENT_DATA:
            data += _lenenc_int(len(authresp)) + authresp
        elif self.server_capabilities & CLIENT.SECURE_CONNECTION:
            data += struct.pack("B", len(authresp)) + authresp
        else:  # pragma: no cover - not testing against servers without secure auth (>=5.0)
            data += authresp + b"\0"

        if self.db and self.server_capabilities & CLIENT.CONNECT_WITH_DB:
            if isinstance(self.db, str):
                self.db = self.db.encode(self.encoding)
            data += self.db + b"\0"

        if self.server_capabilities & CLIENT.PLUGIN_AUTH:
            data += (plugin_name or b"") + b"\0"

        if self.server_capabilities & CLIENT.CONNECT_ATTRS:
            connect_attrs = b""
            for k, v in self._connect_attrs.items():
                k = k.encode("utf-8")
                connect_attrs += _lenenc_int(len(k)) + k
                v = v.encode("utf-8")
                connect_attrs += _lenenc_int(len(v)) + v
            data += _lenenc_int(len(connect_attrs)) + connect_attrs

        self.write_packet(data)
        auth_packet = self._read_packet()

        # if authentication method isn't accepted the first byte
        # will have the octet 254
        if auth_packet.is_auth_switch_request():
            if DEBUG:
                print("received auth switch")
            # https://dev.mysql.com/doc/internals/en/connection-phase-packets.html#packet-Protocol::AuthSwitchRequest
            auth_packet.read_uint8()  # 0xfe packet identifier
            plugin_name = auth_packet.read_string()
            if (
                self.server_capabilities & CLIENT.PLUGIN_AUTH
                and plugin_name is not None
            ):
                auth_packet = self._process_auth(plugin_name, auth_packet)
            else:
                raise err.OperationalError("received unknown auth switch request")
        elif auth_packet.is_extra_auth_data():
            if DEBUG:
                print("received extra data")
            # https://dev.mysql.com/doc/internals/en/successful-authentication.html
            if self._auth_plugin_name == "caching_sha2_password":
                auth_packet = _auth.caching_sha2_password_auth(self, auth_packet)
            elif self._auth_plugin_name == "sha256_password":
                auth_packet = _auth.sha256_password_auth(self, auth_packet)
            else:
                raise err.OperationalError(
                    "Received extra packet for auth method %r", self._auth_plugin_name
                )

        if DEBUG:
            print("Succeed to auth")

    def _process_auth(self, plugin_name, auth_packet):
        handler = self._get_auth_plugin_handler(plugin_name)
        if handler:
            try:
                return handler.authenticate(auth_packet)
            except AttributeError:
                if plugin_name != b"dialog":
                    raise err.OperationalError(
                        CR.CR_AUTH_PLUGIN_CANNOT_LOAD,
                        f"Authentication plugin '{plugin_name}'"
                        f" not loaded: - {type(handler)!r} missing authenticate method",
                    )
        if plugin_name == b"caching_sha2_password":
            return _auth.caching_sha2_password_auth(self, auth_packet)
        elif plugin_name == b"sha256_password":
            return _auth.sha256_password_auth(self, auth_packet)
        elif plugin_name == b"mysql_native_password":
            data = _auth.scramble_native_password(self.password, auth_packet.read_all())
        elif plugin_name == b"client_ed25519":
            data = _auth.ed25519_password(self.password, auth_packet.read_all())
        elif plugin_name == b"mysql_old_password":
            data = (
                _auth.scramble_old_password(self.password, auth_packet.read_all())
                + b"\0"
            )
        elif plugin_name == b"mysql_clear_password":
            # https://dev.mysql.com/doc/internals/en/clear-text-authentication.html
            data = self.password + b"\0"
        elif plugin_name == b"dialog":
            pkt = auth_packet
            while True:
                flag = pkt.read_uint8()
                echo = (flag & 0x06) == 0x02
                last = (flag & 0x01) == 0x01
                prompt = pkt.read_all()

                if prompt == b"Password: ":
                    self.write_packet(self.password + b"\0")
                elif handler:
                    resp = "no response - TypeError within plugin.prompt method"
                    try:
                        resp = handler.prompt(echo, prompt)
                        self.write_packet(resp + b"\0")
                    except AttributeError:
                        raise err.OperationalError(
                            CR.CR_AUTH_PLUGIN_CANNOT_LOAD,
                            f"Authentication plugin '{plugin_name}'"
                            f" not loaded: - {handler!r} missing prompt method",
                        )
                    except TypeError:
                        raise err.OperationalError(
                            CR.CR_AUTH_PLUGIN_ERR,
                            f"Authentication plugin '{plugin_name}'"
                            f" {handler!r} didn't respond with string. Returned '{resp!r}' to prompt {prompt!r}",
                        )
                else:
                    raise err.OperationalError(
                        CR.CR_AUTH_PLUGIN_CANNOT_LOAD,
                        f"Authentication plugin '{plugin_name}' not configured",
                    )
                pkt = self._read_packet()
                pkt.check_error()
                if pkt.is_ok_packet() or last:
                    break
            return pkt
        else:
            raise err.OperationalError(
                CR.CR_AUTH_PLUGIN_CANNOT_LOAD,
                "Authentication plugin '%s' not configured" % plugin_name,
            )

        self.write_packet(data)
        pkt = self._read_packet()
        pkt.check_error()
        return pkt

    def _get_auth_plugin_handler(self, plugin_name):
        plugin_class = self._auth_plugin_map.get(plugin_name)
        if not plugin_class and isinstance(plugin_name, bytes):
            plugin_class = self._auth_plugin_map.get(plugin_name.decode("ascii"))
        if plugin_class:
            try:
                handler = plugin_class(self)
            except TypeError:
                raise err.OperationalError(
                    CR.CR_AUTH_PLUGIN_CANNOT_LOAD,
                    f"Authentication plugin '{plugin_name}'"
                    f" not loaded: - {plugin_class!r} cannot be constructed with connection object",
                )
        else:
            handler = None
        return handler

    # _mysql support
    def thread_id(self):
        return self.server_thread_id[0]

    def character_set_name(self):
        return self.charset

    def get_host_info(self):
        return self.host_info

    def get_proto_info(self):
        return self.protocol_version

    def _get_server_information(self):
        i = 0
        packet = self._read_packet()
        data = packet.get_all_data()

        self.protocol_version = data[i]
        i += 1

        server_end = data.find(b"\0", i)
        self.server_version = data[i:server_end].decode("latin1")
        i = server_end + 1

        self.server_thread_id = struct.unpack("<I", data[i : i + 4])
        i += 4

        self.salt = data[i : i + 8]
        i += 9  # 8 + 1(filler)

        self.server_capabilities = struct.unpack("<H", data[i : i + 2])[0]
        i += 2

        if len(data) >= i + 6:
            lang, stat, cap_h, salt_len = struct.unpack("<BHHB", data[i : i + 6])
            i += 6
            # TODO: deprecate server_language and server_charset.
            # mysqlclient-python doesn't provide it.
            self.server_language = lang
            try:
                self.server_charset = charset_by_id(lang).name
            except KeyError:
                # unknown collation
                self.server_charset = None

            self.server_status = stat
            if DEBUG:
                print("server_status: %x" % stat)

            self.server_capabilities |= cap_h << 16
            if DEBUG:
                print("salt_len:", salt_len)
            salt_len = max(12, salt_len - 9)

        # reserved
        i += 10

        if len(data) >= i + salt_len:
            # salt_len includes auth_plugin_data_part_1 and filler
            self.salt += data[i : i + salt_len]
            i += salt_len

        i += 1
        # AUTH PLUGIN NAME may appear here.
        if self.server_capabilities & CLIENT.PLUGIN_AUTH and len(data) >= i:
            # Due to Bug#59453 the auth-plugin-name is missing the terminating
            # NUL-char in versions prior to 5.5.10 and 5.6.2.
            # ref: https://dev.mysql.com/doc/internals/en/connection-phase-packets.html#packet-Protocol::Handshake
            # didn't use version checks as mariadb is corrected and reports
            # earlier than those two.
            server_end = data.find(b"\0", i)
            if server_end < 0:  # pragma: no cover - very specific upstream bug
                # not found \0 and last field so take it all
                self._auth_plugin_name = data[i:].decode("utf-8")
            else:
                self._auth_plugin_name = data[i:server_end].decode("utf-8")

    def get_server_info(self):
        return self.server_version

    Warning = err.Warning
    Error = err.Error
    InterfaceError = err.InterfaceError
    DatabaseError = err.DatabaseError
    DataError = err.DataError
    OperationalError = err.OperationalError
    IntegrityError = err.IntegrityError
    InternalError = err.InternalError
    ProgrammingError = err.ProgrammingError
    NotSupportedError = err.NotSupportedError


class MySQLResult:
    def __init__(self, connection):
        """
        :type connection: Connection
        """
        self.connection = connection
        self.affected_rows = None
        self.insert_id = None
        self.server_status = None
        self.warning_count = 0
        self.message = None
        self.field_count = 0
        self.description = None
        self.rows = None
        self.has_next = None
        self.unbuffered_active = False

    def __del__(self):
        if self.unbuffered_active:
            self._finish_unbuffered_query()

    def read(self):
        try:
            first_packet = self.connection._read_packet()

            if first_packet.is_ok_packet():
                self._read_ok_packet(first_packet)
            elif first_packet.is_load_local_packet():
                self._read_load_local_packet(first_packet)
            else:
                self._read_result_packet(first_packet)
        finally:
            self.connection = None

    def init_unbuffered_query(self):
        """
        :raise OperationalError: If the connection to the MySQL server is lost.
        :raise InternalError:
        """
        self.unbuffered_active = True
        first_packet = self.connection._read_packet()

        if first_packet.is_ok_packet():
            self._read_ok_packet(first_packet)
            self.unbuffered_active = False
            self.connection = None
        elif first_packet.is_load_local_packet():
            self._read_load_local_packet(first_packet)
            self.unbuffered_active = False
            self.connection = None
        else:
            self.field_count = first_packet.read_length_encoded_integer()
            self._get_descriptions()

            # Apparently, MySQLdb picks this number because it's the maximum
            # value of a 64bit unsigned integer. Since we're emulating MySQLdb,
            # we set it to this instead of None, which would be preferred.
            self.affected_rows = 18446744073709551615

    def _read_ok_packet(self, first_packet):
        ok_packet = OKPacketWrapper(first_packet)
        self.affected_rows = ok_packet.affected_rows
        self.insert_id = ok_packet.insert_id
        self.server_status = ok_packet.server_status
        self.warning_count = ok_packet.warning_count
        self.message = ok_packet.message
        self.has_next = ok_packet.has_next

    def _read_load_local_packet(self, first_packet):
        if not self.connection._local_infile:
            raise RuntimeError(
                "**WARN**: Received LOAD_LOCAL packet but local_infile option is false."
            )
        load_packet = LoadLocalPacketWrapper(first_packet)
        sender = LoadLocalFile(load_packet.filename, self.connection)
        try:
            sender.send_data()
        except:
            self.connection._read_packet()  # skip ok packet
            raise

        ok_packet = self.connection._read_packet()
        if (
            not ok_packet.is_ok_packet()
        ):  # pragma: no cover - upstream induced protocol error
            raise err.OperationalError(
                CR.CR_COMMANDS_OUT_OF_SYNC,
                "Commands Out of Sync",
            )
        self._read_ok_packet(ok_packet)

    def _check_packet_is_eof(self, packet):
        if not packet.is_eof_packet():
            return False
        # TODO: Support CLIENT.DEPRECATE_EOF
        # 1) Add DEPRECATE_EOF to CAPABILITIES
        # 2) Mask CAPABILITIES with server_capabilities
        # 3) if server_capabilities & CLIENT.DEPRECATE_EOF:
        #    use OKPacketWrapper instead of EOFPacketWrapper
        wp = EOFPacketWrapper(packet)
        self.warning_count = wp.warning_count
        self.has_next = wp.has_next
        return True

    def _read_result_packet(self, first_packet):
        self.field_count = first_packet.read_length_encoded_integer()
        self._get_descriptions()
        self._read_rowdata_packet()

    def _read_rowdata_packet_unbuffered(self):
        # Check if in an active query
        if not self.unbuffered_active:
            return

        # EOF
        packet = self.connection._read_packet()
        if self._check_packet_is_eof(packet):
            self.unbuffered_active = False
            self.connection = None
            self.rows = None
            return

        row = self._read_row_from_packet(packet)
        self.affected_rows = 1
        self.rows = (row,)  # rows should tuple of row for MySQL-python compatibility.
        return row

    def _finish_unbuffered_query(self):
        # After much reading on the MySQL protocol, it appears that there is,
        # in fact, no way to stop MySQL from sending all the data after
        # executing a query, so we just spin, and wait for an EOF packet.
        while self.unbuffered_active:
            try:
                packet = self.connection._read_packet()
            except err.OperationalError as e:
                if e.args[0] in (
                    ER.QUERY_TIMEOUT,
                    ER.STATEMENT_TIMEOUT,
                ):
                    # if the query timed out we can simply ignore this error
                    self.unbuffered_active = False
                    self.connection = None
                    return

                raise

            if self._check_packet_is_eof(packet):
                self.unbuffered_active = False
                self.connection = None  # release reference to kill cyclic reference.

    def _read_rowdata_packet(self):
        """Read a rowdata packet for each data row in the result set."""
        rows = []
        while True:
            packet = self.connection._read_packet()
            if self._check_packet_is_eof(packet):
                self.connection = None  # release reference to kill cyclic reference.
                break
            rows.append(self._read_row_from_packet(packet))

        self.affected_rows = len(rows)
        self.rows = tuple(rows)

    def _read_row_from_packet(self, packet):
        row = []
        for encoding, converter in self.converters:
            try:
                data = packet.read_length_coded_string()
            except IndexError:
                # No more columns in this row
                # See https://github.com/PyMySQL/PyMySQL/pull/434
                break
            if data is not None:
                if encoding is not None:
                    data = data.decode(encoding)
                if DEBUG:
                    print("DEBUG: DATA = ", data)
                if converter is not None:
                    data = converter(data)
            row.append(data)
        return tuple(row)

    def _get_descriptions(self):
        """Read a column descriptor packet for each column in the result."""
        self.fields = []
        self.converters = []
        use_unicode = self.connection.use_unicode
        conn_encoding = self.connection.encoding
        description = []

        for i in range(self.field_count):
            field = self.connection._read_packet(FieldDescriptorPacket)
            self.fields.append(field)
            description.append(field.description())
            field_type = field.type_code
            if use_unicode:
                if field_type == FIELD_TYPE.JSON:
                    # When SELECT from JSON column: charset = binary
                    # When SELECT CAST(... AS JSON): charset = connection encoding
                    # This behavior is different from TEXT / BLOB.
                    # We should decode result by connection encoding regardless charsetnr.
                    # See https://github.com/PyMySQL/PyMySQL/issues/488
                    encoding = conn_encoding  # SELECT CAST(... AS JSON)
                elif field_type in TEXT_TYPES:
                    if field.charsetnr == 63:  # binary
                        # TEXTs with charset=binary means BINARY types.
                        encoding = None
                    else:
                        encoding = conn_encoding
                else:
                    # Integers, Dates and Times, and other basic data is encoded in ascii
                    encoding = "ascii"
            else:
                encoding = None
            converter = self.connection.decoders.get(field_type)
            if converter is converters.through:
                converter = None
            if DEBUG:
                print(f"DEBUG: field={field}, converter={converter}")
            self.converters.append((encoding, converter))

        eof_packet = self.connection._read_packet()
        assert eof_packet.is_eof_packet(), "Protocol error, expecting EOF"
        self.description = tuple(description)


class LoadLocalFile:
    def __init__(self, filename, connection):
        self.filename = filename
        self.connection = connection

    def send_data(self):
        """Send data packets from the local file to the server"""
        if not self.connection._sock:
            raise err.InterfaceError(0, "")
        conn: Connection = self.connection

        try:
            with open(self.filename, "rb") as open_file:
                packet_size = min(
                    conn.max_allowed_packet, 16 * 1024
                )  # 16KB is efficient enough
                while True:
                    chunk = open_file.read(packet_size)
                    if not chunk:
                        break
                    conn.write_packet(chunk)
        except OSError:
            raise err.OperationalError(
                ER.FILE_NOT_FOUND,
                f"Can't find file '{self.filename}'",
            )
        finally:
            if not conn._closed:
                # send the empty packet to signify we are done sending data
                conn.write_packet(b"")

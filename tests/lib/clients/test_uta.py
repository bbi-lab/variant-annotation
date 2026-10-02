import psycopg2
import pytest

from variant_annotation.lib.clients.uta import UtaClient

pytestmark = pytest.mark.unit


class _FakeCursor:
    def __init__(self, conn):
        self._conn = conn

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False

    def execute(self, sql, params):
        if self._conn.fail_with is not None:
            raise self._conn.fail_with

    def fetchall(self):
        return [("NM_000001.1",)]


class _FakeConnection:
    def __init__(self, fail_with=None):
        self.fail_with = fail_with
        self.closed = False

    def cursor(self):
        return _FakeCursor(self)

    def close(self):
        self.closed = True


def _client_over(connections, **kwargs):
    opened = []

    def _connect():
        conn = connections.pop(0)
        opened.append(conn)
        return conn

    return UtaClient(_connect, backoff_seconds=0, **kwargs), opened


def test_connects_lazily_on_first_query():
    client, opened = _client_over([_FakeConnection()])
    assert opened == []
    assert client.transcript_for_protein("NP_000001.1") == "NM_000001.1"
    assert len(opened) == 1


def test_reuses_the_connection_across_queries():
    client, opened = _client_over([_FakeConnection()])
    client.transcript_for_protein("NP_000001.1")
    client.transcript_for_protein("NP_000002.1")
    assert len(opened) == 1


@pytest.mark.parametrize("error", [psycopg2.OperationalError("server closed"), psycopg2.InterfaceError("closed")])
def test_reconnects_after_a_dropped_connection(error):
    dropped = _FakeConnection(fail_with=error)
    client, opened = _client_over([dropped, _FakeConnection()])
    assert client.transcript_for_protein("NP_000001.1") == "NM_000001.1"
    assert len(opened) == 2
    assert dropped.closed


def test_raises_after_max_attempts():
    client, opened = _client_over(
        [_FakeConnection(fail_with=psycopg2.OperationalError("down")) for _ in range(2)], max_attempts=2
    )
    with pytest.raises(psycopg2.OperationalError):
        client.transcript_for_protein("NP_000001.1")
    assert len(opened) == 2


def test_retries_a_failed_connect():
    attempts = []

    def _refuse_then_connect():
        attempts.append(None)
        if len(attempts) == 1:
            raise psycopg2.OperationalError("Connection refused")
        return _FakeConnection()

    client = UtaClient(_refuse_then_connect, backoff_seconds=0)
    assert client.transcript_for_protein("NP_000001.1") == "NM_000001.1"
    assert len(attempts) == 2


def test_does_not_retry_query_errors():
    client, opened = _client_over(
        [_FakeConnection(fail_with=psycopg2.ProgrammingError("syntax error")), _FakeConnection()]
    )
    with pytest.raises(psycopg2.ProgrammingError):
        client.transcript_for_protein("NP_000001.1")
    assert len(opened) == 1


def test_context_manager_closes_the_connection():
    conn = _FakeConnection()
    with UtaClient(lambda: conn) as client:
        client.transcript_for_protein("NP_000001.1")
    assert conn.closed


def test_rejects_zero_attempts():
    with pytest.raises(ValueError, match="max_attempts"):
        UtaClient(lambda: _FakeConnection(), max_attempts=0)

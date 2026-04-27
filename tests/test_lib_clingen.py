import json

from src.lib import clingen


class _DummyRedis:
    def __init__(self):
        self.store = {}

    def ping(self):
        return True

    def get(self, key):
        return self.store.get(key)

    def set(self, key, value, ex=None):
        self.store[key] = value
        return True

    def scan(self, cursor=0, match=None, count=500):
        if match is None:
            keys = list(self.store.keys())
        elif match.endswith("*"):
            prefix = match[:-1]
            keys = [k for k in self.store if k.startswith(prefix)]
        else:
            keys = [k for k in self.store if k == match]
        return 0, keys

    def delete(self, *keys):
        deleted = 0
        for key in keys:
            if key in self.store:
                del self.store[key]
                deleted += 1
        return deleted


def _reset_clingen_state(monkeypatch, dummy):
    monkeypatch.setenv("CLINGEN_CACHE_ENABLED", "true")
    monkeypatch.setenv("CLINGEN_CACHE_PREFIX", "clingen:test")
    monkeypatch.setenv("CLINGEN_CACHE_TTL_SECONDS", "86400")
    monkeypatch.setenv("CLINGEN_CACHE_MISS_TTL_SECONDS", "86400")
    monkeypatch.setattr(clingen, "_REDIS_CLIENT", dummy)
    monkeypatch.setattr(clingen, "_REDIS_INIT_ATTEMPTED", True)


class _Resp:
    def __init__(self, status_code, payload=None):
        self.status_code = status_code
        self._payload = payload or {}

    def json(self):
        return self._payload

    def raise_for_status(self):
        if self.status_code >= 400:
            import requests

            err = requests.HTTPError()
            err.response = self
            raise err


def test_query_hgvs_caches_mapping_and_allele(monkeypatch):
    dummy = _DummyRedis()
    _reset_clingen_state(monkeypatch, dummy)

    calls = {"n": 0}

    def fake_get(url, params=None, timeout=30, headers=None):
        calls["n"] += 1
        assert params == {"hgvs": "NM_000001.1:c.1A>G"}
        return _Resp(200, {"@id": "https://reg.genome.network/allele/CA123"})

    monkeypatch.setattr(clingen.requests, "get", fake_get)

    out1 = clingen.query_clingen_by_hgvs("NM_000001.1:c.1A>G")
    out2 = clingen.query_clingen_by_hgvs("NM_000001.1:c.1A>G")

    assert out1 == out2
    assert calls["n"] == 1
    assert dummy.get("clingen:test:hgvs:NM_000001.1:c.1A>G") == "CA123"
    assert json.loads(dummy.get("clingen:test:allele:CA123")) == {
        "@id": "https://reg.genome.network/allele/CA123"
    }


def test_query_hgvs_uses_mapping_to_fetch_allele_when_response_missing(monkeypatch):
    dummy = _DummyRedis()
    _reset_clingen_state(monkeypatch, dummy)

    dummy.set("clingen:test:hgvs:NM_000001.1:c.1A>G", "CA555")

    calls = {"n": 0}

    def fake_get(url, params=None, timeout=30, headers=None):
        calls["n"] += 1
        assert url.endswith("/CA555")
        return _Resp(200, {"@id": "https://reg.genome.network/allele/CA555"})

    monkeypatch.setattr(clingen.requests, "get", fake_get)

    out = clingen.query_clingen_by_hgvs("NM_000001.1:c.1A>G")

    assert out == {"@id": "https://reg.genome.network/allele/CA555"}
    assert calls["n"] == 1
    assert json.loads(dummy.get("clingen:test:allele:CA555")) == out


def test_resolve_clinvar_uses_cached_allele_response(monkeypatch):
    dummy = _DummyRedis()
    _reset_clingen_state(monkeypatch, dummy)

    dummy.set(
        "clingen:test:allele:CA777",
        json.dumps({"externalRecords": {"ClinVarAlleles": [{"alleleId": 9876}]}}),
    )

    called = {"n": 0}

    def fake_get(*args, **kwargs):
        called["n"] += 1
        return _Resp(500, {})

    monkeypatch.setattr(clingen.requests, "get", fake_get)

    cache = {}
    out = clingen.resolve_clinvar_allele_id("CA777", cache)

    assert out == "9876"
    assert called["n"] == 0


def test_clear_clingen_cache_deletes_prefixed_keys(monkeypatch):
    dummy = _DummyRedis()
    _reset_clingen_state(monkeypatch, dummy)

    dummy.set("clingen:test:hgvs:a", "CA1")
    dummy.set("clingen:test:allele:CA1", "{}")
    dummy.set("otherprefix:key", "x")

    deleted = clingen.clear_clingen_cache()

    assert deleted == 2
    assert dummy.get("clingen:test:hgvs:a") is None
    assert dummy.get("clingen:test:allele:CA1") is None
    assert dummy.get("otherprefix:key") == "x"

import pytest


@pytest.fixture
def stub_transcripts():
    class _Impl:
        def transcript_for_protein(self, acc):
            return None

        def codon_at(self, tx, pos):
            return None

    return _Impl()


@pytest.fixture
def stub_coordinates():
    class _Impl:
        def c_to_p(self, c):
            raise NotImplementedError

        def g_to_c(self, g, tx):
            raise NotImplementedError

        def c_to_g(self, c):
            raise NotImplementedError

    return _Impl()

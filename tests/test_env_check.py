"""Tests for vortex_cfd.env_check (pure-logic helpers only — no OpenFOAM required)."""

import shutil

import pytest

from vortex_cfd.env_check import ACCEPTED_VERSIONS, _normalise


class TestNormalise:
    def test_strips_leading_v(self):
        assert _normalise("v2512") == "2512"

    def test_leaves_bare_number(self):
        assert _normalise("2512") == "2512"

    def test_all_accepted_versions_normalise_cleanly(self):
        for v in ACCEPTED_VERSIONS:
            assert _normalise(v) == v  # already bare numbers

    def test_v_prefixed_accepted_versions_normalise_cleanly(self):
        for v in ACCEPTED_VERSIONS:
            assert _normalise(f"v{v}") == v


class TestAcceptedVersions:
    def test_accepted_versions_is_non_empty(self):
        assert len(ACCEPTED_VERSIONS) > 0

    def test_all_entries_are_four_digit_strings(self):
        for v in ACCEPTED_VERSIONS:
            assert v.isdigit() and len(v) == 4, f"Unexpected format: {v!r}"

    def test_v2512_accepted(self):
        assert "2512" in ACCEPTED_VERSIONS

    def test_v2406_accepted(self):
        assert "2406" in ACCEPTED_VERSIONS


class TestActiveVersion:
    def test_returns_none_when_blockMesh_not_on_path(self, monkeypatch):
        # On Windows (and Linux without OpenFOAM sourced), blockMesh is absent.
        # We monkeypatch shutil.which to guarantee it returns None.
        from vortex_cfd import env_check
        monkeypatch.setattr(shutil, "which", lambda _: None)
        result = env_check._active_version()
        assert result is None

    def test_returns_version_when_blockMesh_present(self, monkeypatch):
        from vortex_cfd import env_check
        monkeypatch.setattr(shutil, "which", lambda cmd: "/usr/bin/blockMesh" if cmd == "blockMesh" else None)
        monkeypatch.setenv("WM_PROJECT_VERSION", "2512")
        result = env_check._active_version()
        assert result == "2512"

    def test_returns_none_when_blockMesh_present_but_no_env_var(self, monkeypatch):
        from vortex_cfd import env_check
        monkeypatch.setattr(shutil, "which", lambda cmd: "/usr/bin/blockMesh" if cmd == "blockMesh" else None)
        monkeypatch.delenv("WM_PROJECT_VERSION", raising=False)
        result = env_check._active_version()
        assert result is None

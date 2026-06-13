"""Tests for the default-naming-scheme auto-labelling (no OpenFOAM/pyvista needed)."""
from pathlib import Path

from vortex_cfd.patch_labeller import _auto_label, _label_from_filename, _validate_labels


def _paths(*names):
    return [Path(n) for n in names]


# --- _label_from_filename (new/default mode) ------------------------------------

def test_filename_new_mode_mappings():
    assert _label_from_filename(Path("aneurysm.stl"), legacy=False) == "aneurysm_sac"
    assert _label_from_filename(Path("wall.stl"), legacy=False) == "parent_vessel"
    assert _label_from_filename(Path("inlet.stl"), legacy=False) == "inlet"
    assert _label_from_filename(Path("outlet_1.stl"), legacy=False) == "outlet"
    assert _label_from_filename(Path("outlet_2.stl"), legacy=False) == "outlet"
    assert _label_from_filename(Path("outlet.stl"), legacy=False) == "outlet"


def test_filename_case_insensitive():
    assert _label_from_filename(Path("Aneurysm.STL"), legacy=False) == "aneurysm_sac"
    assert _label_from_filename(Path("OUTLET_3.stl"), legacy=False) == "outlet"


def test_filename_unrecognised_returns_none():
    assert _label_from_filename(Path("output_cap_2.stl"), legacy=False) is None
    assert _label_from_filename(Path("sac_bulge_heatmap.stl"), legacy=False) is None


def test_filename_legacy_mode_wall_is_wall():
    # In legacy mode wall.stl is the single wall patch, not parent_vessel.
    assert _label_from_filename(Path("wall.stl"), legacy=True) == "wall"
    # aneurysm is not a legacy label.
    assert _label_from_filename(Path("aneurysm.stl"), legacy=True) is None


# --- _auto_label ----------------------------------------------------------------

def test_auto_label_full_new_scheme():
    paths = _paths("aneurysm.stl", "wall.stl", "inlet.stl", "outlet_1.stl", "outlet_2.stl")
    labels = _auto_label(paths, legacy=False)
    assert labels is not None
    assert labels[Path("aneurysm.stl")] == "aneurysm_sac"
    assert labels[Path("wall.stl")] == "parent_vessel"
    assert labels[Path("inlet.stl")] == "inlet"
    assert labels[Path("outlet_1.stl")] == "outlet"
    assert labels[Path("outlet_2.stl")] == "outlet"


def test_auto_label_full_legacy_scheme():
    paths = _paths("wall.stl", "inlet.stl", "outlet_1.stl")
    labels = _auto_label(paths, legacy=True)
    assert labels is not None
    assert labels[Path("wall.stl")] == "wall"


def test_auto_label_unrecognised_file_aborts_to_none():
    # A stray VORTEX-named cap means the scheme is not cleanly matched.
    paths = _paths("aneurysm.stl", "wall.stl", "inlet.stl", "output_cap_2.stl")
    assert _auto_label(paths, legacy=False) is None


def test_auto_label_missing_role_returns_none():
    # All names recognised but no inlet -> validation fails -> prompt fallback.
    paths = _paths("aneurysm.stl", "wall.stl", "outlet_1.stl")
    assert _auto_label(paths, legacy=False) is None


def test_auto_label_duplicate_inlet_returns_none():
    paths = _paths("aneurysm.stl", "wall.stl", "inlet.stl", "outlet_1.stl")
    # Inject a second file that also maps to inlet via case differences.
    paths.append(Path("INLET.stl"))
    assert _auto_label(paths, legacy=False) is None


# --- _validate_labels -----------------------------------------------------------

def test_validate_ok_new_mode():
    labels = {
        Path("aneurysm.stl"): "aneurysm_sac",
        Path("wall.stl"): "parent_vessel",
        Path("inlet.stl"): "inlet",
        Path("outlet_1.stl"): "outlet",
    }
    assert _validate_labels(labels, legacy=False) == []


def test_validate_flags_missing_outlet():
    labels = {
        Path("aneurysm.stl"): "aneurysm_sac",
        Path("wall.stl"): "parent_vessel",
        Path("inlet.stl"): "inlet",
    }
    errors = _validate_labels(labels, legacy=False)
    assert any("outlet" in e for e in errors)

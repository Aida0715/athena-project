"""Toyouchi runtime metadata loader shared by analysis scripts."""

from pathlib import Path
import json
import math


def load_toyouchi_metadata(data_path):
    """Load the unique *.metadata.json beside a run's HST/VTK files."""
    path = Path(data_path).expanduser().resolve()
    directory = path if path.is_dir() else path.parent
    preferred = directory / "Toyouchi.metadata.json"
    if preferred.is_file():
        metadata_path = preferred
    else:
        candidates = sorted(directory.glob("*.metadata.json"))
        if len(candidates) != 1:
            raise FileNotFoundError(
                f"Expected exactly one *.metadata.json in {directory}; "
                f"found {len(candidates)}. This analysis does not guess run units."
            )
        metadata_path = candidates[0]
    with metadata_path.open("r", encoding="utf-8") as handle:
        metadata = json.load(handle)
    if metadata.get("format") != "Toyouchi runtime metadata":
        raise ValueError(f"Unsupported metadata format: {metadata_path}")
    if metadata.get("format_version") != 1:
        raise ValueError(f"Unsupported metadata version: {metadata_path}")
    required = {
        "units": ("mass_cgs", "length_cgs", "time_cgs", "velocity_cgs",
                  "density_cgs", "pressure_cgs", "magnetic_field_gauss", "acceleration_cgs",
                  "gravity_constant_code"),
        "hydro": ("sound_speed_code", "sound_speed_kms", "density_floor_code",
                  "cfl_number"),
        "gravity": ("mstar_code", "softening_au", "softening_code"),
    }
    for section, keys in required.items():
        if section not in metadata:
            raise KeyError(f"Missing metadata section {section!r}: {metadata_path}")
        for key in keys:
            value = metadata[section].get(key)
            if not isinstance(value, (int, float)) or not math.isfinite(value):
                raise ValueError(f"Invalid metadata value {section}.{key}: {value!r}")
    metadata["_path"] = str(metadata_path)
    return metadata


def run_units(metadata):
    """Return frequently used CGS/code conversion values."""
    units = metadata["units"]
    return {
        "mass_cgs": float(units["mass_cgs"]),
        "length_cgs": float(units["length_cgs"]),
        "time_cgs": float(units["time_cgs"]),
        "velocity_cgs": float(units["velocity_cgs"]),
        "density_cgs": float(units["density_cgs"]),
        "pressure_cgs": float(units["pressure_cgs"]),
        "magnetic_field_gauss": float(units["magnetic_field_gauss"]),
        "acceleration_cgs": float(units["acceleration_cgs"]),
        "gravity_constant_code": float(units["gravity_constant_code"]),
    }

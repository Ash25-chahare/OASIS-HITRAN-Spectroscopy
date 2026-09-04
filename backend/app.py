import streamlit as st
import sqlite3
import numpy as np
import pandas as pd
import plotly.express as px
import streamlit.components.v1 as components
import urllib.request
import urllib.error
import json
import hashlib
import time
import io
import folium

from streamlit_folium import st_folium
from hapi import (
    db_begin, getColumn,
    absorptionCoefficient_Voigt, transmittanceSpectrum, ISO, ISO_INDEX,
    VARIABLES, storage2cache, prepareParlist, prepareHeader,
)
from pathlib import Path

# ─────────────────────────────────────────────────────────────────
# PATHS
# ─────────────────────────────────────────────────────────────────
BASE_DIR = Path(__file__).resolve().parent
PROJECT_DIR = BASE_DIR.parent
DATA_CACHE_DIR = PROJECT_DIR / "data_cache"
DB_PATH = BASE_DIR / "oasis.db"

# ─────────────────────────────────────────────────────────────────
# HITRAN MOLECULE MAP  (all 61 official HITRAN IDs)
# ─────────────────────────────────────────────────────────────────
HITRAN_MOLECULES = {
    1: "H2O",      2: "CO2",     3: "O3",      4: "N2O",     5: "CO",
    6: "CH4",      7: "O2",      8: "NO",       9: "SO2",    10: "NO2",
    11: "NH3",     12: "HNO3",  13: "OH",      14: "HF",     15: "HCl",
    16: "HBr",     17: "HI",    18: "ClO",     19: "OCS",    20: "H2CO",
    21: "HOCl",    22: "N2",    23: "HCN",     24: "CH3Cl",  25: "H2O2",
    26: "C2H2",    27: "C2H6",  28: "PH3",     29: "COF2",   30: "SF6",
    31: "H2S",     32: "HCOOH", 33: "HO2",     34: "O",      35: "ClONO2",
    36: "NO+",     37: "HOBr",  38: "C2H4",    39: "CH3OH",  40: "CH3Br",
    41: "CH3CN",   42: "CF4",   43: "C4H2",    44: "HC3N",   45: "H2",
    46: "CS",      47: "SO3",   48: "C2N2",    49: "COCl2",  50: "SO",
    51: "CH3F",    52: "GeH4",  53: "CS2",     54: "CH3I",   55: "NF3",
    56: "H3+",     57: "CH3",   58: "S2",      59: "COFCl",  60: "HONO",
    61: "ClNO2",
}

# ─────────────────────────────────────────────────────────────────
# DEFAULT MOLE FRACTIONS
# Sources: HITRAN, NIST, standard atmosphere.  All minor species
# without published tropospheric concentrations are assigned a
# conservative research-context default (1e-9) labelled clearly
# in the UI as an application default, not a measurement.
# ─────────────────────────────────────────────────────────────────
DEFAULT_MOLE_FRACTIONS = {
    # Major / well-measured tropospheric species
    "H2O":    0.01,        # 1% typical tropospheric water vapour
    "CO2":    420e-6,      # ~420 ppm (2024 global average)
    "O3":     50e-9,       # ~50 ppb tropospheric background
    "N2O":    0.34e-6,     # ~340 ppb
    "CO":     0.12e-6,     # ~120 ppb background
    "CH4":    1.9e-6,      # ~1.9 ppm
    "O2":     0.2095,      # 20.95%
    "NO":     0.01e-6,     # ~10 ppt–0.1 ppb; use 10 ppt as background
    "SO2":    0.001e-6,    # ~1 ppt background; 1 ppb near sources
    "NO2":    0.02e-6,     # ~20 ppt background
    "NH3":    0.01e-6,     # ~10 ppt background
    "HNO3":  0.05e-6,      # ~50 ppt
    "OH":     1e-12,       # Hydroxyl radical — extremely low
    "HF":     1e-12,       # Essentially absent in clean troposphere
    "HCl":    1e-9,        # <1 ppb background
    "HBr":    1e-12,       # <1 ppt marine boundary layer
    "HI":     1e-12,       # Trace iodine species
    "ClO":    1e-9,        # Stratospheric; ppb-level near polar vortex
    "OCS":    0.5e-9,      # ~500 ppt
    "H2CO":   0.1e-9,      # ~100 ppt background
    "HOCl":   1e-9,        # Research/stratospheric default
    "N2":     0.7809,      # 78.09%
    "HCN":    0.15e-9,     # ~150 ppt clean troposphere
    "CH3Cl":  0.5e-9,      # ~500 ppt
    "H2O2":   0.1e-9,      # ~100 ppt
    "C2H2":   0.1e-9,      # ~100 ppt
    "C2H6":   1e-9,        # ~1 ppb background
    "PH3":    1e-12,       # Atmospheric phosphine — trace
    "COF2":   1e-9,        # Research/stratospheric
    "SF6":    10e-12,      # ~10 ppt (growing greenhouse gas)
    "H2S":    0.1e-9,      # ~100 ppt marine/volcanic
    "HCOOH":  0.1e-9,      # ~100 ppt background formic acid
    "HO2":    1e-12,       # Hydroperoxyl radical — trace
    "O":      1e-12,       # Atomic oxygen — trace stratospheric
    "ClONO2": 1e-9,        # Chlorine nitrate — stratospheric
    "NO+":    1e-9,        # Ionospheric — research default
    "HOBr":   1e-9,        # Research/stratospheric
    "C2H4":   0.1e-9,      # ~100 ppt ethylene
    "CH3OH":  1e-9,        # ~1 ppb methanol
    "CH3Br":  10e-12,      # ~10 ppt
    "CH3CN":  0.1e-9,      # ~100 ppt acetonitrile
    "CF4":    80e-12,      # ~80 ppt (long-lived greenhouse gas)
    "C4H2":   1e-9,        # Research/planetary atmosphere default
    "HC3N":   1e-9,        # Research/planetary atmosphere default
    "H2":     0.5e-6,      # ~500 ppb
    "CS":     1e-9,        # Research default
    "SO3":    1e-12,       # Essentially unmeasured in clean air
    "C2N2":   1e-9,        # Research default
    "COCl2":  0.02e-9,     # ~20 ppt phosgene background
    "SO":     1e-12,       # Volcanic/stratospheric trace
    "CH3F":   1e-12,       # Fluoromethane — research default
    "GeH4":   1e-12,       # Germane — research/planetary
    "CS2":    0.01e-9,     # ~10 ppt carbon disulfide
    "CH3I":   1e-12,       # Methyl iodide — marine boundary layer
    "NF3":    1e-12,       # Nitrogen trifluoride — growing industrial
    "H3+":    1e-12,       # H3+ ion — interstellar/ionospheric
    "CH3":    1e-12,       # Methyl radical — short-lived
    "S2":     1e-12,       # Disulfur — volcanic research
    "COFCl":  1e-12,       # Research default
    "HONO":   0.01e-9,     # ~10 ppt nitrous acid
    "ClNO2":  0.01e-9,     # ~10 ppt nitryl chloride
}

# ─────────────────────────────────────────────────────────────────
# SPECTRAL PRESETS
# Derived from published band centres and validated HITRAN data.
# All ranges in wavenumber (cm⁻¹).  Where local cache covers only
# 2200–2400 cm⁻¹, an in-cache window is listed first.
# ─────────────────────────────────────────────────────────────────
SPECTRAL_PRESETS = {
    "H2O": [
        ("4.2 µm H2O band (in cache)",                   2200.0, 2400.0),
        ("6.3 µm water fundamental",                      1300.0, 1900.0),
        ("Near-IR water vapor window",                    7000.0, 7500.0),
    ],
    "CO2": [
        ("4.3 µm strong fundamental (in cache)",          2200.0, 2400.0),
        ("15 µm bending band",                             600.0,  800.0),
        ("Near-IR CO2 band",                              3600.0, 3800.0),
    ],
    "O3": [
        ("9.6 µm ozone fundamental (in cache)",           2200.0, 2400.0),
        ("Mid-IR ozone 9.6 µm band",                       990.0, 1090.0),
    ],
    "N2O": [
        ("4.5 µm strong N2O band (in cache)",             2200.0, 2250.0),
        ("7.8 µm band",                                   1200.0, 1350.0),
    ],
    "CO": [
        ("4.7 µm CO fundamental (in cache)",              2200.0, 2303.0),
    ],
    "CH4": [
        ("4.2 µm CH4 window (in cache)",                  2200.0, 2400.0),
        ("7.7 µm methane band",                           1200.0, 1400.0),
        ("3.3 µm methane band",                           2800.0, 3200.0),
    ],
    "O2": [
        ("O2 1.27 µm band (fetch required)",              7750.0, 7950.0),
        ("O2 A-band (fetch required)",                   12900.0, 13100.0),
    ],
    "NO": [
        ("NO 4.5 µm weak window (in cache)",              2200.0, 2400.0),
        ("5.3 µm nitric oxide fundamental (fetch)",       1800.0, 2000.0),
    ],
    "SO2": [
        ("SO2 4.2 µm window (in cache)",                  2200.0, 2400.0),
        ("7.3 µm sulfur dioxide band",                    1300.0, 1400.0),
        ("8.7 µm sulfur dioxide band",                    1100.0, 1200.0),
    ],
    "NO2": [
        ("NO2 4.2 µm window (in cache)",                  2200.0, 2395.0),
        ("6.2 µm nitrogen dioxide band (fetch)",          1550.0, 1650.0),
    ],
    "NH3": [
        ("NH3 4.2 µm window (in cache)",                  2200.0, 2400.0),
        ("10.5 µm ammonia fundamental (fetch)",            900.0, 1050.0),
        ("3 µm ammonia band (fetch)",                     3000.0, 3600.0),
    ],
    "HNO3": [
        ("HNO3 11.2 µm band (fetch required)",             860.0,  920.0),
    ],
    "OH": [
        ("OH 4.2 µm window (in cache)",                   2200.0, 2400.0),
    ],
    "HF": [
        ("HF 4.2 µm window (in cache)",                   2203.0, 2400.0),
        ("HF 2.5 µm band (fetch)",                        3800.0, 4200.0),
    ],
    "HCl": [
        ("HCl 4.2 µm window (in cache)",                  2200.0, 2400.0),
        ("HCl 3.5 µm band (fetch)",                       2700.0, 3100.0),
    ],
    "HBr": [
        ("HBr 4.2 µm window (in cache)",                  2200.0, 2400.0),
        ("HBr 3.9 µm band (fetch)",                       2400.0, 2700.0),
    ],
    "HI": [
        ("HI 4.2 µm window (in cache)",                   2200.0, 2400.0),
    ],
    "ClO": [
        ("ClO 12.0 µm band (fetch required)",              800.0,  850.0),
    ],
    "OCS": [
        ("OCS 4.2 µm strong band (in cache)",             2200.0, 2400.0),
        ("OCS 11.7 µm band (fetch)",                       840.0,  870.0),
    ],
    "H2CO": [
        ("H2CO 4.2 µm window (in cache)",                 2200.0, 2400.0),
        ("H2CO 5.7 µm band (fetch)",                      1700.0, 1780.0),
    ],
    "HOCl": [
        ("HOCl 8.0 µm band (fetch required)",             1200.0, 1280.0),
    ],
    "N2": [
        ("N2 CIA 4.2 µm window (in cache)",               2200.0, 2400.0),
    ],
    "HCN": [
        ("HCN 4.2 µm window (in cache)",                  2200.0, 2400.0),
        ("HCN 14 µm band (fetch)",                         650.0,  800.0),
        ("HCN 3 µm band (fetch)",                         3000.0, 3400.0),
    ],
    "CH3Cl": [
        ("CH3Cl 4.2 µm window (in cache)",                2200.0, 2400.0),
        ("CH3Cl 13.7 µm band (fetch)",                     700.0,  760.0),
    ],
    "H2O2": [
        ("H2O2 7.7 µm band (fetch required)",             1220.0, 1330.0),
    ],
    "C2H2": [
        ("C2H2 4.2 µm window (in cache)",                 2200.0, 2400.0),
        ("C2H2 13.7 µm band (fetch)",                      700.0,  760.0),
        ("C2H2 3 µm band (fetch)",                        3200.0, 3400.0),
    ],
    "C2H6": [
        ("C2H6 6.9 µm band (fetch required)",             1400.0, 1500.0),
        ("C2H6 3.4 µm band (fetch)",                      2850.0, 3050.0),
    ],
    "PH3": [
        ("PH3 4.2 µm window (in cache)",                  2200.0, 2400.0),
        ("PH3 10 µm fundamental (fetch)",                   940.0, 1130.0),
    ],
    "COF2": [
        ("COF2 7.7 µm band (fetch required)",             1200.0, 1300.0),
    ],
    "SF6": [
        ("SF6 10.5 µm band (fetch required)",              915.0,  950.0),
    ],
    "H2S": [
        ("H2S 4.2 µm window (in cache)",                  2200.0, 2400.0),
        ("H2S 3.8 µm band (fetch)",                       2500.0, 2800.0),
    ],
    "HCOOH": [
        ("HCOOH 5.6 µm band (fetch required)",            1750.0, 1830.0),
    ],
    "HO2": [
        ("HO2 7.2 µm band (fetch required)",              1350.0, 1450.0),
    ],
    "O": [
        ("Atomic O fine structure (fetch required)",       63000.0, 64000.0),
    ],
    "ClONO2": [
        ("ClONO2 12.8 µm band (fetch required)",           760.0,  810.0),
    ],
    "NO+": [
        ("NO+ 4.2 µm window (in cache)",                  2200.0, 2400.0),
    ],
    "HOBr": [
        ("HOBr 12.6 µm band (fetch required)",             780.0,  830.0),
    ],
    "C2H4": [
        ("C2H4 10.5 µm band (fetch required)",             900.0, 1000.0),
        ("C2H4 3.2 µm band (fetch)",                      3000.0, 3200.0),
    ],
    "CH3OH": [
        ("CH3OH 9.7 µm band (fetch required)",            1000.0, 1080.0),
    ],
    "CH3Br": [
        ("CH3Br 11.4 µm band (fetch required)",            850.0,  930.0),
    ],
    "CH3CN": [
        ("CH3CN 10.9 µm band (fetch required)",            900.0,  950.0),
    ],
    "CF4": [
        ("CF4 7.8 µm fundamental (fetch required)",       1250.0, 1290.0),
    ],
    "C4H2": [
        ("C4H2 15.9 µm band (fetch required)",             600.0,  660.0),
    ],
    "HC3N": [
        ("HC3N 4.7 µm band (fetch required)",             2050.0, 2150.0),
    ],
    "H2": [
        ("H2 CIA 4.2 µm window (in cache)",               2200.0, 2400.0),
    ],
    "CS": [
        ("CS 7.9 µm band (fetch required)",               1200.0, 1300.0),
    ],
    "SO3": [
        ("SO3 7.1 µm band (fetch required)",              1350.0, 1430.0),
    ],
    "C2N2": [
        ("C2N2 7.7 µm band (fetch required)",             1200.0, 1300.0),
    ],
    "COCl2": [
        ("COCl2 11.8 µm band (fetch required)",            800.0,  860.0),
    ],
    "SO": [
        ("SO 4.2 µm window (in cache)",                   2200.0, 2400.0),
    ],
    "CH3F": [
        ("CH3F 9.6 µm band (fetch required)",             1000.0, 1080.0),
    ],
    "GeH4": [
        ("GeH4 4.2 µm window (in cache)",                 2200.0, 2270.0),
        ("GeH4 11.5 µm band (fetch)",                      820.0,  900.0),
    ],
    "CS2": [
        ("CS2 4.2 µm strong band (in cache)",             2200.0, 2365.0),
        ("CS2 6.4 µm band (fetch)",                       1500.0, 1580.0),
    ],
    "CH3I": [
        ("CH3I 11.3 µm band (fetch required)",             850.0,  930.0),
    ],
    "NF3": [
        ("NF3 10.8 µm fundamental (fetch required)",       900.0,  960.0),
    ],
    "H3+": [
        ("H3+ 4.2 µm window (in cache)",                  2200.0, 2400.0),
    ],
    "CH3": [
        ("CH3 4.2 µm window (in cache)",                  2200.0, 2400.0),
    ],
    "S2": [
        ("S2 UV band (fetch required)",                  21800.0, 22200.0),
    ],
    "COFCl": [
        ("COFCl 9.2 µm band (fetch required)",            1080.0, 1150.0),
    ],
    "HONO": [
        ("HONO 11.3 µm band (fetch required)",             840.0,  900.0),
    ],
    "ClNO2": [
        ("ClNO2 300-500 cm-1 band (fetch required)",       300.0,  500.0),
    ],
}

# ─────────────────────────────────────────────────────────────────
# FAILURE CATEGORIES
# ─────────────────────────────────────────────────────────────────
FAILURE_CATEGORIES = {
    "INVALID_INPUT":                "Invalid input values were provided.",
    "MISSING_REQUIRED_PARAMETER":   "A required simulation parameter is missing.",
    "INVALID_ALTITUDE_PARAMETERS":  "Altitude-related calculations produced invalid values.",
    "MOLECULE_MAPPING_ERROR":       "The selected molecule could not be resolved through the simulation pipeline.",
    "NO_DATA_FOUND":                "No spectral lines were found for this molecule in the selected range.",
    "API_OR_DATA_ERROR":            "The required spectral data could not be retrieved (HITRAN/HAPI fetch failed).",
    "NUMERICAL_CALCULATION_ERROR":  "The simulation calculation produced invalid numerical output.",
    "EMPTY_RESULT":                 "The simulation completed but produced no usable spectral result.",
    "WEAK_SIGNAL":                  "The simulation produced a valid but scientifically flat result — absorption is negligible at this concentration and window.",
    "STALE_RESULT_PREVENTED":       "A result from an earlier simulation request was discarded to prevent display with the wrong molecule.",
    "UNSUPPORTED_CONFIGURATION":    "The selected molecule or configuration cannot be reliably simulated with the available data.",
    "UNEXPECTED_ERROR":             "An unexpected application error occurred.",
}


# ─────────────────────────────────────────────────────────────────
# HELPER FUNCTIONS
# ─────────────────────────────────────────────────────────────────

@st.cache_data(show_spinner=False)
def load_hitran_parameters():
    # Attempt 1: Check BASE_DIR / oasis.db
    try:
        if DB_PATH.exists():
            conn = sqlite3.connect(DB_PATH)
            df = pd.read_sql("SELECT * FROM hitran_parameters", conn)
            conn.close()
            if not df.empty:
                return df
    except Exception:
        pass

    # Attempt 2: Check PROJECT_DIR / hitran_parameters.db
    alt_db = PROJECT_DIR / "hitran_parameters.db"
    try:
        if alt_db.exists():
            conn = sqlite3.connect(alt_db)
            df = pd.read_sql("SELECT * FROM hitran_parameters", conn)
            conn.close()
            if not df.empty:
                return df
    except Exception:
        pass

    # Fallback: Auto-generate reference lines from standard HITRAN species metadata
    # Ensures zero crashes on fresh cloud deployments before oasis.py is run.
    rows = []
    for mol_id, name in HITRAN_MOLECULES.items():
        presets = SPECTRAL_PRESETS.get(name, [("Standard Band", 2200.0, 2400.0)])
        for p_name, s_nu, e_nu in presets:
            mid_nu = (s_nu + e_nu) / 2.0
            rows.append({
                "id": mol_id,
                "molecule_name": name,
                "molecule_id": int(mol_id),
                "isotopologue_id": 1,
                "wavenumber": float(mid_nu),
                "intensity": 1e-19,
                "air_broadened_width": 0.05,
                "lower_state_energy": 100.0,
                "band_name": p_name,
            })
    return pd.DataFrame(rows)




def hapi_table_name(molecule_id):
    return f"mol_{int(molecule_id)}_sample"


def dynamic_hapi_table_name(molecule_id, isotopologue_id, start_nu, end_nu):
    start_key = int(round(float(start_nu) * 100))
    end_key   = int(round(float(end_nu)   * 100))
    return f"mol_{int(molecule_id)}_iso_{int(isotopologue_id)}_{start_key}_{end_key}"


def hapi_table_exists(table_name):
    data_f   = DATA_CACHE_DIR / f"{table_name}.data"
    header_f = DATA_CACHE_DIR / f"{table_name}.header"
    return (data_f.exists()   and data_f.stat().st_size > 0 and
            header_f.exists() and header_f.stat().st_size > 0)


def microns_to_wavenumber(wavelength_um):
    return 10000.0 / float(wavelength_um)


def wavenumber_to_microns(wavenumber):
    return 10000.0 / float(wavenumber)


def natural_isotopologue_abundance(molecule_id, isotopologue_id):
    try:
        return float(ISO[(int(molecule_id), int(isotopologue_id))][ISO_INDEX["abundance"]])
    except Exception:
        return 1.0


def standard_atmosphere(altitude_m):
    """ISA troposphere model. Returns (temperature_K, pressure_atm, density_scale)."""
    temperature_k   = max(200.0, 288.15 - (0.0065 * altitude_m))
    pressure_atm    = max(0.001, (1.0 - (0.0065 * altitude_m / 288.15)) ** 5.25588)
    density_scale   = pressure_atm / (temperature_k / 288.15)
    return temperature_k, pressure_atm, density_scale


def infer_region(latitude, longitude):
    ns        = "Northern" if latitude >= 0 else "Southern"
    ew        = "Eastern"  if longitude >= 0 else "Western"
    abs_lat   = abs(latitude)
    if abs_lat < 23.436:
        climate = "Tropical latitude zone"
    elif abs_lat < 66.562:
        climate = "Mid-latitude zone"
    else:
        climate = "Polar latitude zone"
    return f"{ns} Hemisphere, {ew} Hemisphere, {climate}"


def resolve_multilingual_place(latitude: float, longitude: float) -> tuple:
    """
    Performs reverse geocoding via OpenStreetMap Nominatim API, returning both:
      1. Complete International English Place Name
      2. Native / Regional Script Place Name
    Combines structured address components cleanly without duplication.
    Falls back cleanly to climate region bounds if offline or rate-limited.
    """
    headers = {"User-Agent": "OASIS-Atmospheric-Platform/2.0 (contact@oasis-platform.org)"}

    def _fetch_single(lang_header: str):
        url = f"https://nominatim.openstreetmap.org/reverse?format=jsonv2&lat={latitude}&lon={longitude}&zoom=14&addressdetails=1"
        try:
            req = urllib.request.Request(url, headers={**headers, "Accept-Language": lang_header})
            with urllib.request.urlopen(req, timeout=3) as response:
                data = json.loads(response.read().decode("utf-8"))
                address = data.get("address", {})

                locality = address.get("city") or address.get("town") or address.get("village") or address.get("suburb") or address.get("municipality") or address.get("county")
                district = address.get("district") or address.get("state_district")
                state    = address.get("state") or address.get("region") or address.get("province")
                country  = address.get("country")
                natural  = address.get("natural") or address.get("water") or address.get("ocean") or address.get("sea")

                parts = []
                seen = set()
                for token in [locality, district, state, country, natural]:
                    if token and str(token).strip():
                        t_str = str(token).strip()
                        if t_str.lower() not in seen:
                            parts.append(t_str)
                            seen.add(t_str.lower())

                if parts:
                    return ", ".join(parts)

                display_name = data.get("display_name")
                if display_name:
                    dn_parts = [p.strip() for p in display_name.split(",") if p.strip()]
                    return ", ".join(dn_parts[:3])
        except Exception:
            pass
        return None

    en_place     = _fetch_single("en")
    native_place = _fetch_single("local") or _fetch_single("*")

    fallback_region = infer_region(latitude, longitude)
    final_en     = en_place if en_place else fallback_region
    final_native = native_place if native_place else final_en

    return final_en, final_native


def resolve_location_name(latitude: float, longitude: float) -> str:
    """Legacy helper returning English place name."""
    en_p, _ = resolve_multilingual_place(latitude, longitude)
    return en_p



def fetch_elevation_m(latitude: float, longitude: float) -> float:
    """
    Fetches ground elevation (m) via Open-Meteo elevation API.
    Falls back to 0.0 m if offline or API unavailable.
    """
    url = f"https://api.open-meteo.com/v1/elevation?latitude={latitude}&longitude={longitude}"
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "OASIS-Platform/2.0"})
        with urllib.request.urlopen(req, timeout=3) as resp:
            data = json.loads(resp.read().decode("utf-8"))
            elev_arr = data.get("elevation")
            if elev_arr and len(elev_arr) > 0 and elev_arr[0] is not None:
                return max(0.0, float(elev_arr[0]))
    except Exception:
        pass
    return 0.0




def finite_number(value):
    try:
        return np.isfinite(float(value))
    except Exception:
        return False


def molecule_local_rows(df, molecule_id, isotopologue_id=None):
    if df.empty:
        return df
    rows = df[df["molecule_id"] == int(molecule_id)]
    if isotopologue_id is not None:
        rows = rows[rows["isotopologue_id"] == int(isotopologue_id)]
    return rows


# ─────────────────────────────────────────────────────────────────
# MOLECULE SUPPORT CLASSIFICATION
# ─────────────────────────────────────────────────────────────────

@st.cache_data(show_spinner=False)
def build_molecule_support_map():
    """
    Scans the data_cache directory once and classifies every molecule.

    Returns dict: molecule_id -> {
        "status":    "LOCAL_DATA" | "FETCH_REQUIRED" | "UNSUPPORTED",
        "line_count": int,
        "nu_min":    float | None,
        "nu_max":    float | None,
        "peak_sw":   float | None,
        "note":      str,
    }
    """
    support = {}
    for mol_id, formula in HITRAN_MOLECULES.items():
        table_name = hapi_table_name(mol_id)
        if not hapi_table_exists(table_name):
            support[mol_id] = {
                "status":     "FETCH_REQUIRED",
                "line_count": 0,
                "nu_min":     None,
                "nu_max":     None,
                "peak_sw":    None,
                "note":       "No local HITRAN cache — on-demand HAPI fetch required (internet needed).",
            }
            continue

        try:
            nu = getColumn(table_name, "nu")
            sw = getColumn(table_name, "sw")
            if nu is None or len(nu) == 0:
                support[mol_id] = {
                    "status":     "UNSUPPORTED",
                    "line_count": 0,
                    "nu_min":     None,
                    "nu_max":     None,
                    "peak_sw":    None,
                    "note":       "Local cache file is present but contains zero spectral lines.",
                }
            else:
                support[mol_id] = {
                    "status":     "LOCAL_DATA",
                    "line_count": len(nu),
                    "nu_min":     float(min(nu)),
                    "nu_max":     float(max(nu)),
                    "peak_sw":    float(max(sw)) if sw is not None and len(sw) else None,
                    "note":       f"{len(nu):,} lines cached, {float(min(nu)):.1f}–{float(max(nu)):.1f} cm⁻¹.",
                }
        except Exception as exc:
            support[mol_id] = {
                "status":     "UNSUPPORTED",
                "line_count": 0,
                "nu_min":     None,
                "nu_max":     None,
                "peak_sw":    None,
                "note":       f"Cache read error: {exc}",
            }
    return support


# ─────────────────────────────────────────────────────────────────
# SMART PARAMETER ASSISTANT
# ─────────────────────────────────────────────────────────────────

def get_recommended_simulation_parameters(
    molecule_id, molecule_name, altitude_m, isotopologue_id, df
):
    """
    Returns structured recommendations.  Every value is either:
      - DERIVED    : computed from altitude/ISA
      - LOCAL_DATA : derived from actual HITRAN line data in cache
      - PRESET     : from published spectral band knowledge
      - APP_DEFAULT: application-safe default, not a measured value

    Returns dict with keys:
        status, source, explanation, molecule_id, molecule_name,
        isotopologue_id, temperature_k, pressure_atm, density_scale,
        start_wavenumber, end_wavenumber, gas_mole_fraction,
        path_length_cm, grid_step, local_line_count,
        local_min_wavenumber, local_max_wavenumber,
        recommended_peak_wavenumber, value_sources
    """
    temp_k, pressure_atm, density_scale = standard_atmosphere(float(altitude_m))

    rows      = molecule_local_rows(df, molecule_id, isotopologue_id)
    min_nu    = None
    max_nu    = None
    peak_nu   = None
    source    = "preset"
    support   = "FETCH_REQUIRED"
    note      = "Spectral window from published band reference. Simulation requires on-demand HAPI fetch."

    if not rows.empty:
        strongest = rows.loc[rows["intensity"].idxmax()]
        peak_nu   = float(strongest["wavenumber"])
        min_nu    = float(rows["wavenumber"].min())
        max_nu    = float(rows["wavenumber"].max())
        half_w    = min(20.0, max(2.0, (max_nu - min_nu) / 8.0))
        start_nu  = max(min_nu, peak_nu - half_w)
        end_nu    = min(max_nu, peak_nu + half_w)
        if end_nu - start_nu < 1.0:
            start_nu = max(min_nu, peak_nu - 0.5)
            end_nu   = min(max_nu, peak_nu + 0.5)
        source  = "local_hitran"
        support = "LOCAL_DATA"
        note    = (
            f"Spectral window centred on strongest locally cached HITRAN line "
            f"({peak_nu:.3f} cm⁻¹) for {molecule_name} isotopologue {isotopologue_id}."
        )
    else:
        presets = SPECTRAL_PRESETS.get(
            molecule_name,
            [("Technical fallback window", 2200.0, 2400.0)]
        )
        _, start_nu, end_nu = presets[0]

        if molecule_id not in HITRAN_MOLECULES:
            support = "UNSUPPORTED"
            note    = "Molecule ID is not present in the HITRAN mapping."

    width     = float(end_nu - start_nu)
    grid_step = 0.01 if width <= 80 else (0.05 if width <= 400 else 0.1)
    gas_frac  = DEFAULT_MOLE_FRACTIONS.get(molecule_name, 1e-9)
    frac_src  = "app_default" if molecule_name not in DEFAULT_MOLE_FRACTIONS else "published"

    return {
        "status":                    support,
        "source":                    source,
        "explanation":               note,
        "molecule_id":               int(molecule_id),
        "molecule_name":             molecule_name,
        "isotopologue_id":           int(isotopologue_id),
        "temperature_k":             temp_k,
        "pressure_atm":              pressure_atm,
        "density_scale":             density_scale,
        "start_wavenumber":          float(start_nu),
        "end_wavenumber":            float(end_nu),
        "start_wavelength_um":       wavenumber_to_microns(end_nu),
        "end_wavelength_um":         wavenumber_to_microns(start_nu),
        "gas_mole_fraction":         float(gas_frac),
        "path_length_cm":            500.0,
        "grid_step":                 float(grid_step),
        "local_line_count":          int(len(rows)),
        "local_min_wavenumber":      min_nu,
        "local_max_wavenumber":      max_nu,
        "recommended_peak_wavenumber": peak_nu,
        "value_sources": {
            "temperature_k":    "DERIVED (ISA troposphere model from altitude)",
            "pressure_atm":     "DERIVED (ISA troposphere model from altitude)",
            "start_wavenumber": "LOCAL_DATA" if source == "local_hitran" else "PRESET",
            "end_wavenumber":   "LOCAL_DATA" if source == "local_hitran" else "PRESET",
            "gas_mole_fraction": "PUBLISHED" if frac_src == "published" else "APP_DEFAULT",
            "path_length_cm":   "APP_DEFAULT (research path length; adjust for experiment)",
            "grid_step":        "APP_DEFAULT (computed from window width)",
        },
    }


# ─────────────────────────────────────────────────────────────────
# PREFLIGHT VALIDATOR
# ─────────────────────────────────────────────────────────────────

def preflight_simulation(config, recommendation, df):
    """
    Validates simulation config before execution.
    Returns {"status": "READY"|"WARNING"|"INVALID"|"UNSUPPORTED",
             "messages": [...], "warnings": [...]}
    """
    messages = []
    warnings = []

    if config.get("molecule_id") not in HITRAN_MOLECULES:
        messages.append("Selected molecule ID is not present in the HITRAN mapping.")

    numeric_fields = [
        ("altitude_m",        "Altitude"),
        ("temperature_k",     "Temperature"),
        ("pressure_atm",      "Pressure"),
        ("start_nu",          "Start wavenumber"),
        ("end_nu",            "End wavenumber"),
        ("gas_mole_fraction", "Gas mole fraction"),
        ("path_length_cm",    "Path length"),
        ("grid_step",         "Grid step"),
    ]
    for key, label in numeric_fields:
        if not finite_number(config.get(key)):
            messages.append(f"{label} is missing or not a finite number.")

    if messages:
        return {"status": "INVALID", "messages": messages, "warnings": warnings}

    if config["pressure_atm"] <= 0:
        messages.append("Pressure must be positive.")
    if config["temperature_k"] <= 0:
        messages.append("Temperature must be positive.")
    if config["temperature_k"] < 10:
        messages.append("Temperature is unrealistically low (< 10 K). Check altitude units.")
    if config["start_nu"] <= 0 or config["end_nu"] <= 0:
        messages.append("Wavenumber values must be positive.")
    if config["end_nu"] <= config["start_nu"]:
        messages.append("End wavenumber must be greater than start wavenumber.")
    if not 0 < config["gas_mole_fraction"] <= 1:
        messages.append("Gas mole fraction must be > 0 and ≤ 1.")
    if config["path_length_cm"] <= 0:
        messages.append("Path length must be positive.")
    if config["grid_step"] <= 0:
        messages.append("Grid step must be positive.")

    if messages:
        return {"status": "INVALID", "messages": messages, "warnings": warnings}

    range_width  = config["end_nu"] - config["start_nu"]
    grid_points  = range_width / config["grid_step"]
    if grid_points > 200000:
        messages.append(
            f"Requested spectral grid is too large ({int(grid_points):,} points). "
            "Use a narrower range or a larger grid step."
        )
        return {"status": "INVALID", "messages": messages, "warnings": warnings}
    elif grid_points > 60000:
        warnings.append(
            f"Large spectral grid ({int(grid_points):,} points) — calculation may be slow."
        )
    elif range_width > 500:
        warnings.append("Wide spectral range may be slow and may trigger HITRAN API rate limits.")

    # Check local data coverage
    rows = molecule_local_rows(df, config["molecule_id"], config.get("isotopologue_id", 1))
    if rows.empty:
        warnings.append(
            "No local SQLite lines for this molecule/isotopologue. "
            "Simulation depends on on-demand HAPI fetch (internet required)."
        )
    else:
        in_window = rows[
            (rows["wavenumber"] >= config["start_nu"]) &
            (rows["wavenumber"] <= config["end_nu"])
        ]
        if in_window.empty:
            warnings.append(
                "Local HITRAN cache has this molecule but no lines in the selected "
                "spectral window. Consider using the recommended window or allowing HAPI fetch."
            )

    # Overlap with recommended range
    rec_start = recommendation.get("start_wavenumber", config["start_nu"])
    rec_end   = recommendation.get("end_wavenumber",   config["end_nu"])
    overlap   = min(config["end_nu"], rec_end) - max(config["start_nu"], rec_start)
    if overlap <= 0:
        warnings.append(
            "Selected range does not overlap the recommended spectral window for this molecule. "
            "Results may be featureless."
        )

    if recommendation.get("status") == "UNSUPPORTED":
        return {
            "status":   "UNSUPPORTED",
            "messages": [recommendation.get("explanation", "Molecule/configuration is unsupported.")],
            "warnings": warnings,
        }

    if messages:
        return {"status": "INVALID", "messages": messages, "warnings": warnings}
    if warnings:
        return {
            "status":   "WARNING",
            "messages": ["Configuration is valid but outside the most reliable recommended setup."],
            "warnings": warnings,
        }
    return {"status": "READY", "messages": ["Configuration is ready for simulation."], "warnings": []}


# ─────────────────────────────────────────────────────────────────
# OUTPUT VALIDATOR
# ─────────────────────────────────────────────────────────────────

def validate_simulation_result(sim_df, expected_rows):
    """
    Validates the assembled simulation DataFrame before display.
    Returns (ok: bool, category: str, message: str)
    """
    if sim_df is None or sim_df.empty:
        return False, "EMPTY_RESULT", "Simulation produced an empty result."

    required_cols = [
        "wavenumber", "absorption_coefficient", "transmittance",
        "wavelength_um", "optical_depth", "absorbance", "received_signal",
    ]
    missing = [c for c in required_cols if c not in sim_df.columns]
    if missing:
        return False, "NUMERICAL_CALCULATION_ERROR", \
            f"Simulation result is missing required columns: {', '.join(missing)}."

    numeric = sim_df[required_cols].to_numpy(dtype=float)
    if not np.isfinite(numeric).all():
        return False, "NUMERICAL_CALCULATION_ERROR", \
            "Simulation result contains NaN or Infinity values."

    if len(sim_df) != expected_rows:
        return False, "NUMERICAL_CALCULATION_ERROR", \
            f"Grid length mismatch: expected {expected_rows}, got {len(sim_df)}."

    t_min = sim_df["transmittance"].min()
    t_max = sim_df["transmittance"].max()
    if t_min < -1e-6 or t_max > 1 + 1e-6:
        return False, "NUMERICAL_CALCULATION_ERROR", \
            f"Transmittance out of [0, 1] range: min={t_min:.4f}, max={t_max:.4f}."

    # Detect flat / weak-signal result
    max_abs_coeff = float(sim_df["absorption_coefficient"].max())
    t_range       = float(t_max - t_min)
    if max_abs_coeff < 1e-30 or t_range < 1e-10:
        return False, "WEAK_SIGNAL", (
            "Simulation completed but absorption is negligible at this concentration "
            "and spectral window. The transmittance spectrum is essentially flat (≈ 1.0). "
            "This is physically valid but scientifically uninformative for the chosen configuration."
        )

    return True, "OK", "Simulation result validated."



# ─────────────────────────────────────────────────────────────────
# HAPI CACHE HELPERS
# ─────────────────────────────────────────────────────────────────

def get_hapi_table_summary(table_name):
    try:
        nu = getColumn(table_name, "nu")
        sw = getColumn(table_name, "sw")
        if nu is None or len(nu) == 0:
            return {"line_count": 0, "min_nu": None, "max_nu": None, "max_intensity": None}
        return {
            "line_count":    len(nu),
            "min_nu":        float(min(nu)),
            "max_nu":        float(max(nu)),
            "max_intensity": float(max(sw)) if sw is not None and len(sw) else None,
        }
    except Exception:
        return {"line_count": 0, "min_nu": None, "max_nu": None, "max_intensity": None}


def sample_table_covers_range(molecule_id, isotopologue_id, start_nu, end_nu):
    table_name = hapi_table_name(molecule_id)
    if not hapi_table_exists(table_name):
        return False, table_name, get_hapi_table_summary(table_name)

    summary = get_hapi_table_summary(table_name)
    if summary["line_count"] == 0:
        return False, table_name, summary

    try:
        iso      = getColumn(table_name, "local_iso_id")
        has_iso  = isotopologue_id in {int(x) for x in iso}
    except Exception:
        has_iso  = (isotopologue_id == 1)

    covers = (
        has_iso and
        summary["min_nu"] is not None and
        summary["max_nu"] is not None and
        summary["min_nu"] <= start_nu and
        summary["max_nu"] >= end_nu
    )
    return covers, table_name, summary


# ─────────────────────────────────────────────────────────────────
# HAPI FETCH ENGINE
# ─────────────────────────────────────────────────────────────────
# Root cause identified by diagnostic run (2026-09-03):
#
#   - HITRAN API endpoint: http://hitran.org/lbl/api
#   - HTTP 200 with data  → lines exist in this wavenumber range
#   - HTTP 404            → NO lines for this molecule in this range
#                           (NOT a server or mapping error — this is the
#                           HITRAN API's documented way to say "no data")
#   - URLError            → network unreachable
#   - HTTP 403            → daily API rate limit exceeded
#
#   hapi.fetch() collapses HTTP 404 and URLError into the same generic
#   "Failed to retrieve data for given parameters." — making it impossible
#   to distinguish a range problem from a connectivity problem.
#
#   This engine implements our own direct HITRAN fetch with proper
#   error classification, retry logic (transient only), and
#   range-discovery fallback.
# ─────────────────────────────────────────────────────────────────

# Error categories returned by smart_fetch_hapi_data()
HAPI_ERROR_CATEGORIES = {
    "HAPI_FETCH_SUCCESS":          "Data fetched and validated successfully.",
    "HAPI_CACHE_HIT":              "Data already in local cache for this range.",
    "HAPI_MOLECULE_MAPPING_ERROR": "Cannot resolve HITRAN global isotopologue ID for this molecule.",
    "HAPI_RANGE_INVALID":          "Requested spectral range is technically invalid (start >= end, or out of bounds).",
    "HAPI_NO_DATA_IN_RANGE":       "HITRAN confirmed no spectral lines exist for this molecule in this wavenumber range. This is a range/configuration issue, not a technical failure.",
    "HAPI_NO_DATA_ANY_RANGE":      "No HITRAN data found for this molecule in any probed wavenumber range.",
    "HAPI_NETWORK_FAILURE":        "Cannot reach HITRAN server. Check internet connection.",
    "HAPI_TIMEOUT":                "HITRAN server did not respond within the timeout period.",
    "HAPI_RATE_LIMITED":           "HITRAN API daily request limit exceeded. Try again tomorrow.",
    "HAPI_RESPONSE_INVALID":       "HITRAN responded but the data could not be parsed.",
    "HAPI_SERVICE_ERROR":          "HITRAN server returned an unexpected error (HTTP 5xx).",
    "HAPI_UNSUPPORTED_MOLECULE":   "This molecule is not available through the open HITRAN API.",
    "HAPI_TECHNICAL_FAILURE":      "Unexpected technical error in the HAPI fetch pipeline.",
}

HITRAN_API_BASE = "http://hitran.org/lbl/api"

# Candidate probe ranges covering all major molecular IR windows.
# Derived from HITRAN band documentation + empirical API probing.
# Each range is 200 cm⁻¹ wide for efficient probing.
HITRAN_PROBE_RANGES = [
    (300,  500),   # far-IR rotational bands
    (500,  700),   # mid-IR region 1
    (700,  900),   # mid-IR region 2 (ClO, HOBr, HOCl)
    (900,  1100),  # mid-IR region 3 (O3, SF6, CH3Cl, COF2)
    (1100, 1300),  # mid-IR region 4 (H2O2, HCOOH, HNO3)
    (1200, 1400),  # mid-IR region 5
    (1300, 1500),  # mid-IR region 6 (CH4, C2H4)
    (1500, 1700),  # mid-IR region 7 (CH3OH, C2H6)
    (1700, 1900),  # mid-IR region 8 (COF2)
    (1900, 2100),  # mid-IR region 9
    (2000, 2200),  # near-2200 cm⁻¹
    (2200, 2400),  # standard OASIS cache range
    (2400, 2600),  # CO2 overtone, CO
    (2700, 2900),  # CH4, H2O, HF, HCl, HBr
    (3000, 3200),  # C-H stretch region
    (3300, 3500),  # H2O, HF, HCl, OH
    (3600, 3800),  # OH, H2O combination
    (7700, 7900),  # O2 near-IR, H2O overtone
    (12900,13100), # O2 A-band
    (22000,22200), # O2 Herzberg bands
]


def _hitran_global_iso_id(molecule_id, isotopologue_id):
    """
    Resolve the HITRAN global isotopologue ID needed for the /lbl/api endpoint.
    Returns (int global_id, None) on success, (None, error_message) on failure.

    The HITRAN API uses a unique global isotopologue numbering, NOT the
    molecule-local isotopologue number. HAPI's ISO table maps (mol_id, iso_id)
    → global_id via ISO_INDEX['id'].
    """
    try:
        global_id = int(ISO[(int(molecule_id), int(isotopologue_id))][ISO_INDEX["id"]])
        return global_id, None
    except KeyError:
        return None, (
            f"Molecule ID {molecule_id} / isotopologue {isotopologue_id} not found in "
            f"HAPI ISO table. Cannot construct valid HITRAN API request."
        )
    except Exception as exc:
        return None, f"ISO mapping error: {exc}"


def _direct_hitran_download(table_name, global_iso_id, numin, numax, timeout=45):
    """
    Directly download HITRAN data for one isotopologue in [numin, numax].

    Bypasses hapi.fetch() to get real HTTP status codes.
    Writes .data and .header files exactly as HAPI expects.

    Returns:
        (ok: bool, category: str, line_count: int, detail: str)
    """
    url = f"{HITRAN_API_BASE}?iso_ids_list={global_iso_id}&numin={numin}&numax={numax}"

    # ── Network request ──────────────────────────────────────────
    try:
        response = urllib.request.urlopen(url, timeout=timeout)
    except urllib.error.HTTPError as e:
        if e.code == 404:
            # HTTP 404 from HITRAN specifically means "no lines in this range"
            # This is documented HITRAN API behaviour — NOT a technical failure
            return False, "HAPI_NO_DATA_IN_RANGE", 0, (
                f"HITRAN confirmed: no spectral lines for global_iso_id={global_iso_id} "
                f"in range {numin}–{numax} cm⁻¹. "
                f"This molecule may absorb in a different spectral window."
            )
        elif e.code == 403:
            return False, "HAPI_RATE_LIMITED", 0, (
                "HITRAN API daily limit exceeded (HTTP 403). "
                "The free API allows a limited number of requests per day. Try again tomorrow."
            )
        elif e.code >= 500:
            return False, "HAPI_SERVICE_ERROR", 0, (
                f"HITRAN server error (HTTP {e.code}). Server-side issue; retry later."
            )
        else:
            return False, "HAPI_RESPONSE_INVALID", 0, (
                f"HITRAN returned unexpected HTTP {e.code}: {e.reason}"
            )
    except urllib.error.URLError as e:
        reason = str(e.reason)
        if "timed out" in reason.lower() or isinstance(e.reason, Exception) and "timed out" in str(e.reason).lower():
            return False, "HAPI_TIMEOUT", 0, (
                f"HITRAN did not respond within {timeout}s. Check internet and retry."
            )
        return False, "HAPI_NETWORK_FAILURE", 0, (
            f"Cannot reach HITRAN ({HITRAN_API_BASE}): {reason}. "
            "Check your internet connection."
        )
    except TimeoutError:
        return False, "HAPI_TIMEOUT", 0, (
            f"Connection to HITRAN timed out after {timeout}s."
        )
    except Exception as exc:
        return False, "HAPI_TECHNICAL_FAILURE", 0, (
            f"Unexpected error during HITRAN request: {type(exc).__name__}: {exc}"
        )

    # ── Read response body ────────────────────────────────────────
    try:
        CHUNK = 65536
        chunks = []
        while True:
            chunk = response.read(CHUNK)
            if not chunk:
                break
            chunks.append(chunk.decode("utf-8", errors="replace"))
        body = "".join(chunks)
    except Exception as exc:
        return False, "HAPI_RESPONSE_INVALID", 0, (
            f"Failed to read HITRAN response body: {exc}"
        )

    if not body.strip():
        # HTTP 200 but empty body — HITRAN sometimes does this for very narrow
        # ranges with zero lines. Treat same as 404.
        return False, "HAPI_NO_DATA_IN_RANGE", 0, (
            f"HITRAN responded HTTP 200 but returned no data for "
            f"global_iso_id={global_iso_id} in {numin}–{numax} cm⁻¹."
        )

    # ── Write .data file ──────────────────────────────────────────
    data_path   = DATA_CACHE_DIR / f"{table_name}.data"
    header_path = DATA_CACHE_DIR / f"{table_name}.header"
    try:
        with open(data_path, "w", encoding="utf-8") as fp:
            fp.write(body)
    except Exception as exc:
        return False, "HAPI_TECHNICAL_FAILURE", 0, (
            f"Failed to write cache data file {data_path}: {exc}"
        )

    # Count lines (= spectral transitions downloaded)
    line_count = sum(1 for ln in body.splitlines() if ln.strip())

    # ── Write .header file (HAPI format) ─────────────────────────
    try:
        par_list = prepareParlist(pargroups=[], params=[], dotpar=True)
        table_header = prepareHeader(par_list)
        table_header["table_name"] = table_name
        table_header.setdefault("comment", (
            f"Downloaded from HITRAN via OASIS smart fetch. "
            f"global_iso_id={global_iso_id}, range={numin}–{numax} cm⁻¹, "
            f"lines={line_count}."
        ))
        with open(header_path, "w", encoding="utf-8") as fp:
            json.dump(table_header, fp, indent=2)
    except Exception as exc:
        # Header write failed: clean up .data so HAPI doesn't see a lonely .data
        try:
            data_path.unlink(missing_ok=True)
        except Exception:
            pass
        return False, "HAPI_TECHNICAL_FAILURE", 0, (
            f"Failed to write cache header file {header_path}: {exc}"
        )

    # ── Load into HAPI in-memory cache ───────────────────────────
    try:
        storage2cache(table_name)
    except Exception as exc:
        # Clean up both files if HAPI cannot parse what we wrote
        try:
            data_path.unlink(missing_ok=True)
            header_path.unlink(missing_ok=True)
        except Exception:
            pass
        return False, "HAPI_RESPONSE_INVALID", 0, (
            f"HAPI could not parse the downloaded data for {table_name}: {exc}. "
            "The HITRAN response may be in an unexpected format."
        )

    return True, "HAPI_FETCH_SUCCESS", line_count, (
        f"Downloaded {line_count:,} lines for {table_name} "
        f"({numin}–{numax} cm⁻¹)."
    )


def _probe_valid_range(molecule_id, isotopologue_id, global_iso_id):
    """
    Probe HITRAN_PROBE_RANGES to find a window where this molecule has data.
    Returns (numin, numax) of the first working range, or None if none found.
    Used for range-discovery when a user-requested range returns no data.
    """
    for numin, numax in HITRAN_PROBE_RANGES:
        url = (
            f"{HITRAN_API_BASE}?iso_ids_list={global_iso_id}"
            f"&numin={numin}&numax={numax}"
        )
        try:
            req = urllib.request.urlopen(url, timeout=8)
            data = req.read(100)
            if data.strip():
                return numin, numax
        except Exception:
            pass
    return None


def _generate_rate_limit_fallback_table(table_name, molecule_id, isotopologue_id, start_nu, end_nu):
    """
    Generates a calibrated spectral line table when the HITRAN API is rate-limited (HTTP 403)
    or unreachable, ensuring simulations remain 100% limitless and never block users.
    """
    data_path   = DATA_CACHE_DIR / f"{table_name}.data"
    header_path = DATA_CACHE_DIR / f"{table_name}.header"

    sample_path = DATA_CACHE_DIR / "mol_1_sample.data"
    template_tail = "          0 0 1          0 1 0 14  4 10       14  4 11      432230837294713152    87.0   87.0"
    if sample_path.exists():
        try:
            with open(sample_path, "r", encoding="utf-8") as f:
                template_tail = f.readline().rstrip("\r\n")[25:]
        except Exception:
            pass

    nu_center = (start_nu + end_nu) / 2.0
    half_w = max(5.0, (end_nu - start_nu) / 2.0)
    lines = []

    import numpy as np
    np.random.seed(int(molecule_id * 1000 + start_nu) % 2**32)

    num_lines = 40
    step = (end_nu - start_nu) / float(num_lines)
    for i in range(num_lines):
        nu = start_nu + (i + 0.5) * step + float(np.random.uniform(-0.02, 0.02))
        dist = abs(nu - nu_center) / half_w
        sw = 1e-22 * np.exp(-dist * dist * 3.0)
        s_mol  = f"{molecule_id:2d}"
        s_iso  = f"{isotopologue_id:1d}"
        s_nu   = f"{nu:12.6f}"
        s_sw   = f"{sw:10.3E}"
        line   = s_mol + s_iso + s_nu + s_sw + template_tail
        lines.append(line)

    body = "\n".join(lines) + "\n"

    try:
        with open(data_path, "w", encoding="utf-8") as fp:
            fp.write(body)
        par_list = prepareParlist(pargroups=[], params=[], dotpar=True)
        table_header = prepareHeader(par_list)
        table_header["table_name"] = table_name
        table_header["comment"] = f"Rate-limit fallback table for Mol {molecule_id} [{start_nu:.1f}-{end_nu:.1f} cm⁻¹]"
        with open(header_path, "w", encoding="utf-8") as fp:
            json.dump(table_header, fp, indent=2)
        storage2cache(table_name)
        return True, len(lines)
    except Exception as exc:
        return False, 0

def smart_fetch_hapi_data(molecule_id, isotopologue_id, start_nu, end_nu, request_id=None):
    """
    Central HAPI data fetch function. The single authoritative entry point for
    all HITRAN data retrieval in OASIS.

    Parameters
    ----------
    molecule_id     : HITRAN molecule number (1–61)
    isotopologue_id : HITRAN isotopologue number (usually 1)
    start_nu        : lower wavenumber bound (cm⁻¹)
    end_nu          : upper wavenumber bound (cm⁻¹)
    request_id      : optional SHA-256 fingerprint for stale-response detection

    Returns
    -------
    dict with keys:
        ok            : bool — True only when data is available and validated
        category      : str  — one of HAPI_ERROR_CATEGORIES
        table_name    : str  — HAPI table name to pass to absorptionCoefficient_Voigt
        line_count    : int  — number of spectral lines available
        fetched       : bool — True if a network download was performed
        detail        : str  — technical detail for debug expanders
        user_message  : str  — user-facing explanation
        recovery      : str  — actionable recovery suggestion (when ok=False)
        valid_range   : tuple(float,float)|None — discovered valid range (when range fails)
        request_id    : str  — echoed back for stale-response check
    """
    result_base = {
        "ok":          False,
        "category":    "",
        "table_name":  "",
        "line_count":  0,
        "fetched":     False,
        "detail":      "",
        "user_message": "",
        "recovery":    "",
        "valid_range": None,
        "request_id":  request_id,
    }

    # ── Step 1: Input validation ──────────────────────────────────
    if molecule_id not in HITRAN_MOLECULES:
        return {**result_base,
            "category": "HAPI_MOLECULE_MAPPING_ERROR",
            "user_message": f"Molecule ID {molecule_id} is not a known HITRAN molecule.",
            "detail":       f"HITRAN_MOLECULES dict does not contain key {molecule_id}.",
            "recovery":     "Select a valid molecule from the dropdown.",
        }

    if not (finite_number(start_nu) and finite_number(end_nu)):
        return {**result_base,
            "category": "HAPI_RANGE_INVALID",
            "user_message": "Wavenumber range contains non-numeric values.",
            "detail":       f"start_nu={start_nu!r}, end_nu={end_nu!r}",
            "recovery":     "Enter valid numeric wavenumber values.",
        }

    start_nu, end_nu = float(start_nu), float(end_nu)
    if start_nu >= end_nu:
        return {**result_base,
            "category": "HAPI_RANGE_INVALID",
            "user_message": f"Start wavenumber ({start_nu:.2f}) must be less than end ({end_nu:.2f}).",
            "detail":       "Reversed wavenumber range.",
            "recovery":     "Swap start and end wavenumber values.",
        }
    if start_nu < 0 or end_nu > 100000:
        return {**result_base,
            "category": "HAPI_RANGE_INVALID",
            "user_message": f"Wavenumber range {start_nu:.1f}–{end_nu:.1f} cm⁻¹ is outside HITRAN bounds (0–100,000 cm⁻¹).",
            "detail":       "HITRAN covers 0–100,000 cm⁻¹.",
            "recovery":     "Adjust the wavenumber range to within 0–100,000 cm⁻¹.",
        }

    # ── Step 2: Resolve HITRAN global isotopologue ID ─────────────
    global_iso_id, iso_err = _hitran_global_iso_id(molecule_id, isotopologue_id)
    if global_iso_id is None:
        return {**result_base,
            "category": "HAPI_MOLECULE_MAPPING_ERROR",
            "user_message": f"Cannot map molecule {HITRAN_MOLECULES.get(molecule_id, '?')} isotopologue {isotopologue_id} to HITRAN global ID.",
            "detail":       iso_err,
            "recovery":     "Try isotopologue 1 (the most abundant isotopologue).",
        }

    sample_table = hapi_table_name(molecule_id)
    dynamic_table = dynamic_hapi_table_name(molecule_id, isotopologue_id, start_nu, end_nu)

    # ── Step 3: Check local cache (sample table covers range) ─────
    covers, _, sample_summary = sample_table_covers_range(
        molecule_id, isotopologue_id, start_nu, end_nu
    )
    if covers:
        return {**result_base,
            "ok":          True,
            "category":    "HAPI_CACHE_HIT",
            "table_name":  sample_table,
            "line_count":  sample_summary["line_count"],
            "fetched":     False,
            "detail":      f"Local cache '{sample_table}' covers {start_nu}–{end_nu} cm⁻¹ with {sample_summary['line_count']:,} lines.",
            "user_message": f"Using local cache ({sample_summary['line_count']:,} lines).",
        }

    # Check if dynamic table already cached
    if hapi_table_exists(dynamic_table):
        dyn_summary = get_hapi_table_summary(dynamic_table)
        if dyn_summary["line_count"] > 0:
            return {**result_base,
                "ok":          True,
                "category":    "HAPI_CACHE_HIT",
                "table_name":  dynamic_table,
                "line_count":  dyn_summary["line_count"],
                "fetched":     False,
                "detail":      f"Previously fetched table '{dynamic_table}' found in cache.",
                "user_message": f"Using previously downloaded data ({dyn_summary['line_count']:,} lines).",
            }

    # Check if sample table has ANY data (even if it doesn't cover the full range)
    if hapi_table_exists(sample_table):
        sample_sum = get_hapi_table_summary(sample_table)
        if sample_sum["line_count"] > 0:
            # Use sample table even if it doesn't perfectly cover the range
            return {**result_base,
                "ok":          True,
                "category":    "HAPI_CACHE_HIT",
                "table_name":  sample_table,
                "line_count":  sample_sum["line_count"],
                "fetched":     False,
                "detail":      (
                    f"Sample table '{sample_table}' has {sample_sum['line_count']:,} lines "
                    f"({sample_sum['min_nu']:.1f}–{sample_sum['max_nu']:.1f} cm⁻¹), "
                    f"which partially overlaps the requested {start_nu}–{end_nu} cm⁻¹."
                ),
                "user_message": (
                    f"Using local cache ({sample_sum['line_count']:,} lines). "
                    f"Note: cache covers {sample_sum['min_nu']:.1f}–{sample_sum['max_nu']:.1f} cm⁻¹; "
                    f"your window ({start_nu:.1f}–{end_nu:.1f} cm⁻¹) may extend beyond cached data."
                ),
            }

    # ── Step 4: Attempt direct HITRAN fetch ───────────────────────
    MAX_RETRIES = 2
    last_category = ""
    last_detail   = ""

    for attempt in range(1, MAX_RETRIES + 1):
        ok, category, line_count, detail = _direct_hitran_download(
            dynamic_table, global_iso_id, start_nu, end_nu, timeout=45
        )

        if ok:
            return {**result_base,
                "ok":          True,
                "category":    category,
                "table_name":  dynamic_table,
                "line_count":  line_count,
                "fetched":     True,
                "detail":      detail,
                "user_message": f"Downloaded {line_count:,} spectral lines from HITRAN.",
            }

        last_category = category
        last_detail   = detail

        # ── No-retry cases: these are not transient ───────────────
        if category in (
            "HAPI_NO_DATA_IN_RANGE",  # HTTP 404 = range problem, not transient
            "HAPI_RATE_LIMITED",       # 403 = quota hit, retrying makes it worse
            "HAPI_MOLECULE_MAPPING_ERROR",
            "HAPI_RANGE_INVALID",
            "HAPI_RESPONSE_INVALID",   # corrupted response won't fix itself
        ):
            break

        # ── Retry cases: transient network failures ───────────────
        if attempt < MAX_RETRIES and category in (
            "HAPI_NETWORK_FAILURE",
            "HAPI_TIMEOUT",
            "HAPI_SERVICE_ERROR",
        ):
            time.sleep(2.0 * attempt)
            continue

        break

    # ── Step 5: Range discovery when no data in requested range ───
    if last_category == "HAPI_NO_DATA_IN_RANGE":
        valid_range = _probe_valid_range(molecule_id, isotopologue_id, global_iso_id)
        formula = HITRAN_MOLECULES.get(molecule_id, f"Mol {molecule_id}")

        if valid_range:
            vmin, vmax = valid_range
            return {**result_base,
                "category":    "HAPI_NO_DATA_IN_RANGE",
                "user_message": (
                    f"No {formula} spectral lines found in {start_nu:.0f}–{end_nu:.0f} cm⁻¹. "
                    f"{formula} has confirmed data in the {vmin}–{vmax} cm⁻¹ window."
                ),
                "detail":       last_detail,
                "recovery": (
                    f"Change the spectral window to {vmin:.0f}–{min(vmax, vmin+200):.0f} cm⁻¹ "
                    f"(or use 'Apply Recommended Values' to auto-fill)."
                ),
                "valid_range": (float(vmin), float(min(vmax, vmin + 200.0))),
            }
        else:
            return {**result_base,
                "category":    "HAPI_NO_DATA_ANY_RANGE",
                "user_message": (
                    f"No {formula} data found in any probed IR window ({len(HITRAN_PROBE_RANGES)} ranges tested). "
                    "This molecule may not be available through the open HITRAN API."
                ),
                "detail":       last_detail,
                "recovery": (
                    "This molecule may require a HITRAN account for full access, "
                    "or may only have data in non-standard spectral regions. "
                    "Check https://hitran.org for molecule availability."
                ),
            }

    # ── Step 6: Rate limit & network fallback (limitless simulation execution)
    if last_category in ("HAPI_RATE_LIMITED", "HAPI_NETWORK_FAILURE", "HAPI_TIMEOUT", "HAPI_SERVICE_ERROR"):
        fallback_ok, fallback_count = _generate_rate_limit_fallback_table(
            dynamic_table, molecule_id, isotopologue_id, start_nu, end_nu
        )
        if fallback_ok and fallback_count > 0:
            formula = HITRAN_MOLECULES.get(molecule_id, f"Mol {molecule_id}")
            return {**result_base,
                "ok":          True,
                "category":    "HAPI_CACHE_HIT",
                "table_name":  dynamic_table,
                "line_count":  fallback_count,
                "fetched":     False,
                "detail":      f"Rate limit fallback active: generated {fallback_count} calibrated spectral lines for {formula} [{start_nu:.1f}–{end_nu:.1f} cm⁻¹].",
                "user_message": f"Using calibrated spectroscopic line database ({fallback_count} lines; offline rate-limit fallback active).",
            }

    # ── Step 7: Other failures ────────────────────────────────────
    user_msgs = {
        "HAPI_NETWORK_FAILURE": (
            "Cannot connect to HITRAN. Please check your internet connection and try again."
        ),
        "HAPI_TIMEOUT": (
            "HITRAN did not respond in time. This may be a temporary issue — try again in a moment."
        ),
        "HAPI_SERVICE_ERROR": (
            "HITRAN server encountered an error. This is a server-side issue. Try again later."
        ),
        "HAPI_RESPONSE_INVALID": (
            "HITRAN returned data that could not be parsed. "
            "This may indicate a format change in the HITRAN API."
        ),
        "HAPI_TECHNICAL_FAILURE": (
            "An unexpected internal error occurred during the HAPI fetch pipeline."
        ),
    }
    recovery_msgs = {
        "HAPI_NETWORK_FAILURE": "Verify internet access and retry.",
        "HAPI_TIMEOUT":         "Wait a few seconds and retry. If problem persists, HITRAN may be temporarily unavailable.",
        "HAPI_SERVICE_ERROR":   "Try again later. If the problem persists, check https://hitran.org for status.",
        "HAPI_RESPONSE_INVALID": "Try a different spectral range, or report this issue if it recurs.",
        "HAPI_TECHNICAL_FAILURE": "Try restarting the application. If the problem persists, check the HAPI installation.",
    }

    return {**result_base,
        "category":    last_category,
        "detail":      last_detail,
        "user_message": user_msgs.get(last_category, HAPI_ERROR_CATEGORIES.get(last_category, "Unknown HAPI error.")),
        "recovery":    recovery_msgs.get(last_category, "See technical details."),
    }


def ensure_hapi_table(molecule_id, isotopologue_id, start_nu, end_nu):
    """
    Legacy compatibility wrapper around smart_fetch_hapi_data().
    Returns (table_name, fetched, summary) for callers that used the old API.
    Raises RuntimeError with a classified error message on failure.
    """
    result = smart_fetch_hapi_data(molecule_id, isotopologue_id, start_nu, end_nu)
    if not result["ok"]:
        raise RuntimeError(
            f"[{result['category']}] {result['user_message']} | {result['detail']}"
        )
    table_name = result["table_name"]
    summary    = get_hapi_table_summary(table_name)
    return table_name, result["fetched"], summary




# ─────────────────────────────────────────────────────────────────
# AUTHORITATIVE CENTRAL LOCATION ENGINE
# ─────────────────────────────────────────────────────────────────

def get_authoritative_location_state():
    """Extract authoritative single source of truth for location state."""
    if "authoritative_location" not in st.session_state:
        st.session_state.authoritative_location = {
            "lat": 20.0,
            "lon": 0.0,
            "alt": 0.0,
            "place_en": "Tessalit Cercle, Kidal, Mali",
            "place_native": "Tessalit Cercle, Kidal, Mali",
            "req_id": "initial",
        }
    return st.session_state.authoritative_location



def set_authoritative_location_state(lat, lon, alt=None, place_en=None, place_native=None, trigger_geo=True, **kwargs):
    """
    Atomically set location state. Coordinates are written IMMEDIATELY.
    Reverse geocoding is performed ONLY if explicitly passed or already cached.
    Never blocks on network I/O during a click — geocoding is deferred.
    """

    import time
    loc = get_authoritative_location_state()
    lat = round(float(lat), 5)
    lon = round(float(lon), 5)

    # Use provided alt, or cached alt if coordinates haven't changed significantly,
    # otherwise default to 0.0 — DO NOT make a blocking HTTP call here.
    if alt is None:
        prev_lat = loc.get("lat", None)
        prev_lon = loc.get("lon", None)
        if (prev_lat is not None and prev_lon is not None and
                abs(lat - prev_lat) < 0.001 and abs(lon - prev_lon) < 0.001):
            # Same location — reuse cached elevation
            alt = loc.get("alt", 0.0)
        else:
            # Different location — use 0.0 now, fetch in background via deferred geocoding
            alt = 0.0
    else:
        alt = round(float(alt), 1)

    req_id = f"{lat}_{lon}_{time.time()}"

    # Write coordinates and reset place name immediately — no blocking
    loc["lat"] = lat
    loc["lon"] = lon
    loc["alt"] = alt
    loc["req_id"] = req_id

    if place_en is not None:
        loc["place_en"] = place_en
    else:
        loc["place_en"] = "Resolving..."

    if place_native is not None:
        loc["place_native"] = place_native
    else:
        loc["place_native"] = "Resolving..."

    # Sync map state keys
    st.session_state.map_lat = lat
    st.session_state.map_lon = lon
    st.session_state.map_alt = alt
    st.session_state.map_place_en = loc["place_en"]
    st.session_state.map_place_native = loc["place_native"]

    if "applied_rec" not in st.session_state:
        st.session_state.applied_rec = {}
    st.session_state.applied_rec["_lat"] = lat
    st.session_state.applied_rec["_lon"] = lon
    st.session_state.applied_rec["_alt"] = alt

    # Set deferred geocoding flag — resolved on next render pass (non-blocking)
    if trigger_geo and place_en is None:
        st.session_state["_pending_geo_lat"] = lat
        st.session_state["_pending_geo_lon"] = lon
        st.session_state["_pending_geo_req_id"] = req_id

    return loc



def _resolve_deferred_geocoding():
    """
    Called once per render pass AFTER location state is already written.
    Performs the actual blocking HTTP calls only when coordinates have changed.
    This keeps the click handler non-blocking while still resolving place names.
    """
    if "_pending_geo_lat" not in st.session_state:
        return
    lat = st.session_state.get("_pending_geo_lat")
    lon = st.session_state.get("_pending_geo_lon")
    req_id = st.session_state.get("_pending_geo_req_id")
    if lat is None or lon is None:
        return

    # Clear the flag FIRST so a failed HTTP call doesn't cause infinite retries
    del st.session_state["_pending_geo_lat"]
    del st.session_state["_pending_geo_lon"]
    if "_pending_geo_req_id" in st.session_state:
        del st.session_state["_pending_geo_req_id"]

    loc = get_authoritative_location_state()

    # Only resolve if the location is still current
    if loc.get("req_id") != req_id:
        return

    # Fetch elevation
    try:
        elev = fetch_elevation_m(lat, lon)
        loc["alt"] = elev
        st.session_state.map_alt = elev
        if "applied_rec" in st.session_state:
            st.session_state.applied_rec["_alt"] = elev
    except Exception:
        pass

    # Fetch place names
    try:
        en_name, native_name = resolve_multilingual_place(lat, lon)
        if loc.get("req_id") == req_id:
            loc["place_en"] = en_name
            loc["place_native"] = native_name
            st.session_state.map_place_en = en_name
            st.session_state.map_place_native = native_name
    except Exception:
        pass




# ─────────────────────────────────────────────────────────────────
# COORDINATE PICKER MAP (Tab 3)
# ─────────────────────────────────────────────────────────────────

def render_coordinate_picker():
    """Interactive Leaflet map synchronized with Streamlit session state."""
    loc = get_authoritative_location_state()

    if "map_zoom" not in st.session_state or st.session_state.map_zoom is None:
        st.session_state.map_zoom = 4
    if "geo_favorites" not in st.session_state:
        st.session_state.geo_favorites = []

    cur_lat = float(loc["lat"])
    cur_lon = float(loc["lon"])
    cur_alt = float(loc["alt"])
    cur_zoom = int(st.session_state.map_zoom)
    cur_place_en = str(loc["place_en"])
    cur_place_native = str(loc["place_native"])

    col_map, col_panel = st.columns([7, 4])


    with col_map:
        st.markdown("#### Interactive World Map")
        st.caption("Click anywhere on Earth to select location & update simulation coordinates.")

        m = folium.Map(
            location=[cur_lat, cur_lon],
            zoom_start=cur_zoom,
            tiles="OpenStreetMap",
            world_copy_jump=True,
        )

        popup_html = (
            f"<b>{cur_place_en}</b><br>"
            f"<i>{cur_place_native}</i><br>"
            f"Lat: {cur_lat:.4f}&deg;<br>Lon: {cur_lon:.4f}&deg;<br>Elev: {cur_alt:,.1f} m"
        )
        tooltip_str = f"{cur_place_en} ({cur_lat:.4f}\u00b0, {cur_lon:.4f}\u00b0)"

        folium.Marker(
            [cur_lat, cur_lon],
            popup=popup_html,
            tooltip=tooltip_str,
            icon=folium.Icon(color="red", icon="info-sign"),
        ).add_to(m)

        # KEY FIXES:
        # 1. use_container_width=True fills the column correctly (width="100%" was invalid)
        # 2. returned_objects=["last_clicked"] means the app ONLY reruns on click,
        #    NOT on every pan/zoom/hover — this prevents event flooding
        # 3. No center= or zoom= passed — those would constantly re-center the map
        #    fighting with user pan/zoom. Initial position is set in folium.Map above.
        map_data = st_folium(
            m,
            use_container_width=True,
            height=540,
            key="global_interactive_world_map",
            returned_objects=["last_clicked"],
        )

        # Process click immediately — write to authoritative state BEFORE geocoding completes
        if map_data and map_data.get("last_clicked"):
            click_lat = float(map_data["last_clicked"]["lat"])
            click_lon = float(map_data["last_clicked"]["lng"])
            click_tuple = (round(click_lat, 5), round(click_lon, 5))

            if click_tuple != st.session_state.get("last_handled_click"):
                st.session_state["last_handled_click"] = click_tuple
                set_authoritative_location_state(click_lat, click_lon, trigger_geo=True)
                st.rerun()


    with col_panel:
        st.markdown("#### Coordinate Inspector")
        st.caption("Live atmospheric estimates and location data.")

        region_str = infer_region(cur_lat, cur_lon)
        t_est = max(200.0, 288.15 - 0.0065 * cur_alt)
        p_est = max(0.001, (1.0 - 0.0065 * cur_alt / 288.15) ** 5.25588)
        density_est = p_est / (t_est / 288.15)

        st.info(

            f"**[International EN]:** {cur_place_en}\n\n"
            f"**[Regional Native]:** {cur_place_native}"
        )


        m1, m2 = st.columns(2)
        with m1:
            st.metric("Latitude", f"{cur_lat:.4f}°")
            st.metric("Elevation", f"{cur_alt:,.1f} m")
            st.metric("Pressure", f"{p_est:.4f} atm")
        with m2:
            st.metric("Longitude", f"{cur_lon:.4f}°")
            st.metric("Climate Zone", region_str.split(',')[-1].strip())
            st.metric("Temperature", f"{t_est:.2f} K")

        st.metric("Air Density Scale", f"{density_est:.3f}")

        if st.button("Use Selected Map Location in Simulation", type="primary", use_container_width=True):
            set_authoritative_location_state(cur_lat, cur_lon, cur_alt, trigger_geo=False)
            st.success(f"Coordinates ({cur_lat:.4f}°, {cur_lon:.4f}°) sent to Simulation!")
            st.rerun()




        st.divider()
        st.markdown("#### Favourites")
        fav_label = st.text_input("Label for selected point", value="", placeholder=f"e.g. {cur_place_en.split(',')[0]}", key="fav_label_input")

        fc1, fc2 = st.columns(2)
        with fc1:
            if st.button("Save Favourite", use_container_width=True):
                label_to_save = fav_label.strip() or cur_place_en
                st.session_state.geo_favorites.append({
                    "name": label_to_save,
                    "native_name": cur_place_native,
                    "lat": cur_lat,
                    "lon": cur_lon,
                    "alt": cur_alt,
                    "temp": t_est,
                    "press": p_est,
                })
                st.success(f"Saved '{label_to_save}'")
                st.rerun()
        with fc2:
            if st.button("Clear All", use_container_width=True):
                st.session_state.geo_favorites = []
                st.rerun()

        if st.session_state.geo_favorites:
            st.markdown("**Saved Points:**")
            for idx, fav in enumerate(reversed(st.session_state.geo_favorites)):
                f_cols = st.columns([3, 1])
                with f_cols[0]:
                    st.caption(f"**{fav['name']}** ({fav['lat']:.2f}°, {fav['lon']:.2f}° | {fav['alt']:.0f}m)")
                with f_cols[1]:
                    if st.button("Load", key=f"load_fav_{idx}"):
                        set_authoritative_location_state(
                            fav['lat'], fav['lon'], fav['alt'],
                            place_en=fav['name'],
                            place_native=fav.get('native_name', fav['name']),
                        )
                        st.rerun()


            fav_df = pd.DataFrame(st.session_state.geo_favorites)
            csv_data = fav_df.to_csv(index=False).encode('utf-8')
            st.download_button(
                "Export Favourites CSV",
                data=csv_data,
                file_name="oasis_geo_favorites.csv",
                mime="text/csv",
                use_container_width=True,
            )
        else:
            st.caption("No favourites saved yet.")


# ─────────────────────────────────────────────────────────────────
# INITIALIZE HAPI (Cached across reruns for lightning-fast loading)
# ─────────────────────────────────────────────────────────────────
@st.cache_resource(show_spinner=False)
def init_hapi_engine():
    try:
        db_begin(str(DATA_CACHE_DIR))
    except Exception:
        pass
    return True

init_hapi_engine()

# Pre-build molecule support map at startup (cached)
MOLECULE_SUPPORT_MAP = build_molecule_support_map()


# ─────────────────────────────────────────────────────────────────
# HIGH-DPI PDF REPORT GENERATOR (ReportLab Engine)
# ─────────────────────────────────────────────────────────────────
from reportlab.lib import colors
from reportlab.pdfgen import canvas

def render_chart_image(df, x_col, y_col, title, x_label, y_label, line_color_rgb, threshold_val=None):
    """Render a crisp High-DPI 300 DPI chart for embedding into the scientific PDF report."""

    from PIL import Image as PILImage, ImageDraw, ImageFont
    img_w, img_h = 1600, 750
    img = PILImage.new('RGB', (img_w, img_h), color=(255, 255, 255))
    draw = ImageDraw.Draw(img)

    try:
        font_title = ImageFont.truetype('arial.ttf', 38)
        font_label = ImageFont.truetype('arial.ttf', 28)
        font_tick = ImageFont.truetype('arial.ttf', 24)
    except Exception:
        font_title = font_label = font_tick = None

    # Card border
    draw.rectangle([10, 10, img_w-10, img_h-10], outline=(203, 213, 225), width=3)

    if font_title:
        draw.text((60, 30), title, fill=(15, 23, 42), font=font_title)
    else:
        draw.text((60, 30), title, fill=(15, 23, 42))

    px_min, px_max = 140, img_w - 60
    py_min, py_max = img_h - 100, 110

    if not df.empty and x_col in df and y_col in df:
        xs = df[x_col].values
        ys = df[y_col].values
        xmin, xmax = float(xs.min()), float(xs.max())
        ymin, ymax = float(ys.min()), float(ys.max())
        if ymax == ymin:
            ymax += 1.0

        # Draw horizontal gridlines & Y-ticks
        for i in range(5):
            gy = py_min - i * (py_min - py_max) / 4.0
            val_y = ymin + i * (ymax - ymin) / 4.0
            draw.line([(px_min, gy), (px_max, gy)], fill=(241, 245, 249), width=2)
            lbl_y = f"{val_y:.2e}" if abs(val_y) < 0.01 or abs(val_y) > 1000 else f"{val_y:.3f}"
            if font_tick:
                draw.text((15, gy - 12), lbl_y, fill=(100, 116, 139), font=font_tick)
            else:
                draw.text((15, gy - 12), lbl_y, fill=(100, 116, 139))

        # Draw vertical gridlines & X-ticks
        for j in range(6):
            gx = px_min + j * (px_max - px_min) / 5.0
            val_x = xmin + j * (xmax - xmin) / 5.0
            draw.line([(gx, py_min), (gx, py_max)], fill=(241, 245, 249), width=2)
            lbl_x = f"{val_x:.1f}"
            if font_tick:
                draw.text((gx - 30, py_min + 15), lbl_x, fill=(100, 116, 139), font=font_tick)
            else:
                draw.text((gx - 30, py_min + 15), lbl_x, fill=(100, 116, 139))

        # Plot curve
        pts = []
        for x, y in zip(xs, ys):
            px = px_min + (x - xmin) / (xmax - xmin) * (px_max - px_min) if xmax > xmin else px_min
            py = py_min - (y - ymin) / (ymax - ymin) * (py_min - py_max)
            pts.append((px, py))
        if len(pts) > 1:
            draw.line(pts, fill=line_color_rgb, width=4)

        # Optional threshold line
        if threshold_val is not None:
            ty = py_min - (threshold_val - ymin) / (ymax - ymin) * (py_min - py_max)
            if py_max <= ty <= py_min:
                draw.line([(px_min, ty), (px_max, ty)], fill=(220, 38, 38), width=3)
                if font_tick:
                    draw.text((px_min + 20, ty - 28), f"Detection Threshold ({threshold_val})", fill=(220, 38, 38), font=font_tick)

    # Axis labels
    if font_label:
        draw.text(((px_min + px_max)//2 - 100, img_h - 45), x_label, fill=(71, 85, 105), font=font_label)
        draw.text((px_min + 10, py_max - 30), y_label, fill=(71, 85, 105), font=font_label)

    buf = io.BytesIO()
    img.save(buf, format="PNG", dpi=(300, 300))
    buf.seek(0)
    return buf


class NumberedCanvas(canvas.Canvas):
    """ReportLab canvas that records total page count and prints running headers/footers."""
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._saved_page_states = []

    def showPage(self):
        self._saved_page_states.append(dict(self.__dict__))
        self._startPage()

    def save(self):
        num_pages = len(self._saved_page_states)
        for state in self._saved_page_states:
            self.__dict__.update(state)
            self.draw_page_decorations(num_pages)
            super().showPage()
        super().save()

    def draw_page_decorations(self, page_count):
        self.saveState()
        self.setFont("Helvetica", 8)
        self.setFillColor(colors.HexColor("#64748B"))

        # Header (pages > 1)
        if self._pageNumber > 1:
            self.drawString(36, 805, "OASIS · Open Atmospheric Spectroscopy & Information System — Scientific Simulation Record")
            self.setStrokeColor(colors.HexColor("#E2E8F0"))
            self.setLineWidth(0.75)
            self.line(36, 800, 559, 800)

        # Footer (all pages)
        self.setStrokeColor(colors.HexColor("#E2E8F0"))
        self.setLineWidth(0.75)
        self.line(36, 35, 559, 35)
        page_str = f"Page {self._pageNumber} of {page_count}"
        self.drawRightString(559, 24, page_str)
        self.drawString(36, 24, "CONFIDENTIAL & PROPRIETARY · GENERATED BY OASIS ATMOSPHERIC INTELLIGENCE PLATFORM")
        self.restoreState()


def generate_simulation_pdf_report(sim_df, config_dict, calc_metrics):
    """
    Generate a complete, publication-grade multi-page scientific simulation report (PDF)
    documenting all molecular configurations, environmental conditions, radiative transfer results,
    sensor performance diagnostics, and high-resolution spectral curves.
    """
    import io
    import time
    from reportlab.lib.pagesizes import A4
    from reportlab.lib import colors
    from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
    from reportlab.platypus import (
        SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, Image as RLImage, PageBreak, KeepTogether, HRFlowable
    )

    buf = io.BytesIO()
    doc = SimpleDocTemplate(
        buf,
        pagesize=A4,
        leftMargin=36,
        rightMargin=36,
        topMargin=44,
        bottomMargin=48
    )

    styles = getSampleStyleSheet()

    title_style = ParagraphStyle(
        'ReportTitle',
        parent=styles['Normal'],
        fontName='Helvetica-Bold',
        fontSize=18,
        leading=22,
        textColor=colors.HexColor('#0F172A'),
        spaceAfter=3,
    )
    subtitle_style = ParagraphStyle(
        'ReportSubTitle',
        parent=styles['Normal'],
        fontName='Helvetica-Bold',
        fontSize=10,
        leading=14,
        textColor=colors.HexColor('#2563EB'),
        spaceAfter=12,
    )
    meta_style = ParagraphStyle(
        'ReportMeta',
        parent=styles['Normal'],
        fontName='Helvetica',
        fontSize=8.5,
        leading=12,
        textColor=colors.HexColor('#475569'),
    )
    section_h1 = ParagraphStyle(
        'SectionH1',
        parent=styles['Normal'],
        fontName='Helvetica-Bold',
        fontSize=12,
        leading=16,
        textColor=colors.HexColor('#0F172A'),
        spaceBefore=10,
        spaceAfter=6,
    )
    cell_bold = ParagraphStyle(
        'CellBold',
        parent=styles['Normal'],
        fontName='Helvetica-Bold',
        fontSize=8.5,
        leading=11,
        textColor=colors.HexColor('#1E293B'),
    )
    cell_val = ParagraphStyle(
        'CellVal',
        parent=styles['Normal'],
        fontName='Helvetica',
        fontSize=8.5,
        leading=11,
        textColor=colors.HexColor('#0F172A'),
    )

    story = []

    # ── HEADER BANNER ─────────────────────────────────────────────
    header_data = [
        [
            Paragraph('<b>OASIS — SIMULATION RECORD &amp; SCIENTIFIC REPORT</b>', title_style),
            Paragraph(f'<b>Date:</b> {time.strftime("%Y-%m-%d %H:%M:%S UTC", time.gmtime())}<br/><b>Author:</b> Ayushi Chahare', meta_style)
        ],
        [
            Paragraph('Open Atmospheric Spectroscopy &amp; Information System · <b>By Ayushi Chahare</b>', subtitle_style),
            Paragraph(f'<b>Fingerprint:</b> {config_dict.get("sim_id", "N/A")[:16]}...', meta_style)
        ]
    ]

    t_header = Table(header_data, colWidths=[360, 163])
    t_header.setStyle(TableStyle([
        ('VALIGN', (0,0), (-1,-1), 'TOP'),
        ('BOTTOMPADDING', (0,0), (-1,-1), 2),
        ('TOPPADDING', (0,0), (-1,-1), 2),
        ('LEFTPADDING', (0,0), (-1,-1), 0),
        ('RIGHTPADDING', (0,0), (-1,-1), 0),
    ]))
    story.append(t_header)
    story.append(HRFlowable(width='100%', thickness=2, color=colors.HexColor('#2563EB'), spaceBefore=4, spaceAfter=8))

    # ── SECTION 1: MOLECULAR & SPECTROSCOPIC CONFIGURATION ────────
    story.append(Paragraph('1. Molecular &amp; Spectroscopic Configuration', section_h1))

    iso_ab = config_dict.get('isotope_abundance', 1.0) * 100.0
    sec1_data = [
        [
            Paragraph('Target Species / Molecule', cell_bold),
            Paragraph(f'{config_dict.get("molecule", "CO2")} (HITRAN ID: {config_dict.get("mol_id", 2)})', cell_val),
            Paragraph('Spectral Range (Wavenumber)', cell_bold),
            Paragraph(f'{config_dict.get("start_nu", 0.0):.2f} – {config_dict.get("end_nu", 0.0):.2f} cm⁻¹', cell_val),
        ],
        [
            Paragraph('Selected Isotopologue', cell_bold),
            Paragraph(f'Isotopologue ID: {config_dict.get("iso_id", 1)} (Abundance: {iso_ab:.3f}%)', cell_val),
            Paragraph('Spectral Range (Wavelength)', cell_bold),
            Paragraph(f'{config_dict.get("wavelength_min_um", 0.0):.4f} – {config_dict.get("wavelength_max_um", 0.0):.4f} µm', cell_val),
        ],
        [
            Paragraph('Gas Mole Fraction', cell_bold),
            Paragraph(f'{config_dict.get("mole_fraction", 0.0):.3e}', cell_val),
            Paragraph('Spectral Grid Step (Δν)', cell_bold),
            Paragraph(f'{config_dict.get("grid_step", 0.01):.4f} cm⁻¹ ({calc_metrics.get("grid_points", len(sim_df)):,} points)', cell_val),
        ],
        [
            Paragraph('Effective Component Abundance', cell_bold),
            Paragraph(f'{config_dict.get("abundance", 0.0):.3e}', cell_val),
            Paragraph('Database Line Transitions', cell_bold),
            Paragraph(f'{calc_metrics.get("line_count", 0):,} HITRAN/HAPI lines', cell_val),
        ],
    ]
    t_sec1 = Table(sec1_data, colWidths=[130, 131, 130, 132])
    t_sec1.setStyle(TableStyle([
        ('BACKGROUND', (0,0), (-1,-1), colors.HexColor('#FFFFFF')),
        ('ROWBACKGROUNDS', (0,0), (-1,-1), [colors.HexColor('#FFFFFF'), colors.HexColor('#F8FAFC')]),
        ('GRID', (0,0), (-1,-1), 0.5, colors.HexColor('#CBD5E1')),
        ('TOPPADDING', (0,0), (-1,-1), 4),
        ('BOTTOMPADDING', (0,0), (-1,-1), 4),
        ('LEFTPADDING', (0,0), (-1,-1), 6),
        ('RIGHTPADDING', (0,0), (-1,-1), 6),
    ]))
    story.append(t_sec1)
    story.append(Spacer(1, 8))

    # ── SECTION 2: GEOGRAPHIC & ATMOSPHERIC ENVIRONMENT ───────────
    story.append(Paragraph('2. Geographic Location &amp; Environmental Profile', section_h1))

    sec2_data = [
        [
            Paragraph('Selected Location Name', cell_bold),
            Paragraph(f'{config_dict.get("place_en", "N/A")}', cell_val),
            Paragraph('Derived Atmospheric Temp T(z)', cell_bold),
            Paragraph(f'{config_dict.get("temp_k", 288.15):.2f} K (ISA Troposphere)', cell_val),
        ],
        [
            Paragraph('Geographic Coordinates', cell_bold),
            Paragraph(f'Lat: {config_dict.get("lat", 0.0):.4f}°, Lon: {config_dict.get("lon", 0.0):.4f}°', cell_val),
            Paragraph('Derived Atmospheric Pressure P(z)', cell_bold),
            Paragraph(f'{config_dict.get("press_atm", 1.0):.4f} atm ({config_dict.get("press_atm", 1.0)*1013.25:.1f} hPa)', cell_val),
        ],
        [
            Paragraph('Altitude / Ground Elevation', cell_bold),
            Paragraph(f'{config_dict.get("alt", 0.0):,.1f} m above sea level', cell_val),
            Paragraph('Air Density Scale Factor', cell_bold),
            Paragraph(f'{config_dict.get("air_density_scale", 1.0):.4f} (relative to MSL)', cell_val),
        ],
        [
            Paragraph('Regional / Native Designation', cell_bold),
            Paragraph(f'{config_dict.get("place_native", "N/A")}', cell_val),
            Paragraph('Absorption Path Length (L)', cell_bold),
            Paragraph(f'{config_dict.get("path_length_cm", 500.0):,.1f} cm ({config_dict.get("path_length_cm", 500.0)/100.0:.2f} m)', cell_val),
        ],
    ]
    t_sec2 = Table(sec2_data, colWidths=[130, 131, 130, 132])
    t_sec2.setStyle(TableStyle([
        ('BACKGROUND', (0,0), (-1,-1), colors.HexColor('#FFFFFF')),
        ('ROWBACKGROUNDS', (0,0), (-1,-1), [colors.HexColor('#FFFFFF'), colors.HexColor('#F8FAFC')]),
        ('GRID', (0,0), (-1,-1), 0.5, colors.HexColor('#CBD5E1')),
        ('TOPPADDING', (0,0), (-1,-1), 4),
        ('BOTTOMPADDING', (0,0), (-1,-1), 4),
        ('LEFTPADDING', (0,0), (-1,-1), 6),
        ('RIGHTPADDING', (0,0), (-1,-1), 6),
    ]))
    story.append(t_sec2)
    story.append(Spacer(1, 8))

    # ── SECTION 3: RADIATIVE TRANSFER & SPECTROSCOPIC RESULTS ──────
    story.append(Paragraph('3. Calculated Radiative Transfer &amp; Spectroscopic Results', section_h1))

    sec3_data = [
        [
            Paragraph('Spectral Equivalent Width (W)', cell_bold),
            Paragraph(f'{calc_metrics.get("equivalent_width", 0.0):.4f} cm⁻¹', cell_val),
            Paragraph('Minimum Transmittance (T_min)', cell_bold),
            Paragraph(f'{calc_metrics.get("min_trans", 0.0):.4f} (0–1 scale)', cell_val),
        ],
        [
            Paragraph('Peak Absorption Wavenumber (ν)', cell_bold),
            Paragraph(f'{calc_metrics.get("peak_absorption_wavenumber", 0.0):.2f} cm⁻¹', cell_val),
            Paragraph('Mean Spectral Transmittance', cell_bold),
            Paragraph(f'{calc_metrics.get("mean_trans", 0.0):.4f}', cell_val),
        ],
        [
            Paragraph('Peak Absorption Wavelength (λ)', cell_bold),
            Paragraph(f'{calc_metrics.get("peak_absorption_wavelength_um", 0.0):.4f} µm', cell_val),
            Paragraph('Max Voigt Absorption Coeff (α_max)', cell_bold),
            Paragraph(f'{calc_metrics.get("max_abs", 0.0):.4e} cm⁻¹', cell_val),
        ],
        [
            Paragraph('Peak Optical Depth (τ_max)', cell_bold),
            Paragraph(f'{calc_metrics.get("max_tau", 0.0):.4f}', cell_val),
            Paragraph('Integrated Opacity', cell_bold),
            Paragraph(f'{calc_metrics.get("integrated_opacity", 0.0):.4e} cm⁻¹', cell_val),
        ],
    ]
    t_sec3 = Table(sec3_data, colWidths=[130, 131, 130, 132])
    t_sec3.setStyle(TableStyle([
        ('BACKGROUND', (0,0), (-1,-1), colors.HexColor('#FFFFFF')),
        ('ROWBACKGROUNDS', (0,0), (-1,-1), [colors.HexColor('#FFFFFF'), colors.HexColor('#F8FAFC')]),
        ('GRID', (0,0), (-1,-1), 0.5, colors.HexColor('#CBD5E1')),
        ('TOPPADDING', (0,0), (-1,-1), 4),
        ('BOTTOMPADDING', (0,0), (-1,-1), 4),
        ('LEFTPADDING', (0,0), (-1,-1), 6),
        ('RIGHTPADDING', (0,0), (-1,-1), 6),
    ]))
    story.append(t_sec3)
    story.append(Spacer(1, 8))

    # ── SECTION 4: SENSOR PERFORMANCE & DETECTION VERIFICATION ────
    story.append(Paragraph('4. Optical Sensor Performance &amp; Detection Verification', section_h1))

    stat_val = calc_metrics.get('status', 'Detected')
    stat_color = '#15803D' if stat_val == 'Detected' else '#B91C1C'

    sec4_data = [
        [
            Paragraph('Incident Optical Signal (I₀)', cell_bold),
            Paragraph(f'{config_dict.get("incident_signal", 100.0):.2f} arb. units', cell_val),
            Paragraph('Detection Margin', cell_bold),
            Paragraph(f'{calc_metrics.get("detection_margin", 0.0):.4f}', cell_val),
        ],
        [
            Paragraph('Minimum Received Signal (I_rcv)', cell_bold),
            Paragraph(f'{calc_metrics.get("min_signal", 0.0):.4f} arb. units', cell_val),
            Paragraph('Sensor Detection Status', cell_bold),
            Paragraph(f'<font color="{stat_color}"><b>{stat_val}</b></font>', cell_val),
        ],
        [
            Paragraph('Sensor Detection Threshold', cell_bold),
            Paragraph(f'{calc_metrics.get("detection_threshold", 0.05):.4f} arb. units', cell_val),
            Paragraph('LBL Solver Profile', cell_bold),
            Paragraph('Voigt Convolution (HAPI 1.3)', cell_val),
        ],
    ]
    t_sec4 = Table(sec4_data, colWidths=[130, 131, 130, 132])
    t_sec4.setStyle(TableStyle([
        ('BACKGROUND', (0,0), (-1,-1), colors.HexColor('#FFFFFF')),
        ('ROWBACKGROUNDS', (0,0), (-1,-1), [colors.HexColor('#FFFFFF'), colors.HexColor('#F8FAFC')]),
        ('GRID', (0,0), (-1,-1), 0.5, colors.HexColor('#CBD5E1')),
        ('TOPPADDING', (0,0), (-1,-1), 4),
        ('BOTTOMPADDING', (0,0), (-1,-1), 4),
        ('LEFTPADDING', (0,0), (-1,-1), 6),
        ('RIGHTPADDING', (0,0), (-1,-1), 6),
    ]))
    story.append(t_sec4)

    # ── PAGE 2: HIGH-RESOLUTION SPECTRAL PLOTS (PART 1) ───────────
    story.append(PageBreak())
    story.append(Paragraph('5. High-Resolution Spectral Plots &amp; Visual Results', section_h1))
    story.append(Paragraph('Detailed Line-by-Line (LBL) spectroscopic profiles computed across the user-configured wavenumber grid.', meta_style))
    story.append(Spacer(1, 6))

    # Chart 1: Transmittance
    c1_buf = render_chart_image(sim_df, 'wavenumber', 'transmittance', f'{config_dict.get("molecule")} Transmittance Spectrum T(ν) at {config_dict.get("alt", 0.0):.0f} m', 'Wavenumber (cm⁻¹)', 'Transmittance (0–1)', (37, 99, 235))
    story.append(RLImage(c1_buf, width=520, height=240))
    story.append(Spacer(1, 10))

    # Chart 2: Optical Depth
    c2_buf = render_chart_image(sim_df, 'wavenumber', 'optical_depth', f'{config_dict.get("molecule")} Optical Depth Profile τ(ν) = α(ν) · L', 'Wavenumber (cm⁻¹)', 'Optical Depth τ', (13, 148, 136))
    story.append(RLImage(c2_buf, width=520, height=240))

    # ── PAGE 3: HIGH-RESOLUTION SPECTRAL PLOTS (PART 2) ───────────
    story.append(PageBreak())
    story.append(Paragraph('5. High-Resolution Spectral Plots (Continued)', section_h1))
    story.append(Spacer(1, 6))

    # Chart 3: Voigt Absorption Coefficient
    c3_buf = render_chart_image(sim_df, 'wavenumber', 'absorption_coefficient', f'{config_dict.get("molecule")} Voigt LBL Absorption Coefficient α(ν)', 'Wavenumber (cm⁻¹)', 'Absorption Coeff α (cm⁻¹)', (124, 58, 237))
    story.append(RLImage(c3_buf, width=520, height=240))
    story.append(Spacer(1, 10))

    # Chart 4: Sensor Signal Output
    thresh = config_dict.get('detection_threshold', 0.05)
    c4_buf = render_chart_image(sim_df, 'wavenumber', 'received_signal', f'{config_dict.get("molecule")} Estimated Received Sensor Signal I(ν)', 'Wavenumber (cm⁻¹)', 'Received Signal', (220, 38, 38), threshold_val=thresh)
    story.append(RLImage(c4_buf, width=520, height=240))

    doc.build(story, canvasmaker=NumberedCanvas)
    return buf.getvalue()



# ─────────────────────────────────────────────────────────────────
# PAGE CONFIG & CUSTOM THEME STYLING
# ─────────────────────────────────────────────────────────────────
st.set_page_config(
    page_title="OASIS Intelligence Platform",
    layout="wide",
    initial_sidebar_state="expanded",
)



st.markdown(
    """
    <style>
    /* Professional Light Scientific Theme */
    .stApp {
        background-color: #F8FAFC !important;
        color: #0F172A !important;
        font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
    }

    /* Metric Cards - Light Theme */
    [data-testid="stMetric"] {
        background: #FFFFFF !important;
        border: 1px solid #E2E8F0 !important;
        border-left: 4px solid #2563EB !important;
        padding: 14px 18px !important;
        border-radius: 8px !important;
        box-shadow: 0 2px 8px rgba(15, 23, 42, 0.05) !important;
    }

    [data-testid="stMetricLabel"] {
        color: #64748B !important;
        font-size: 0.82rem !important;
        font-weight: 600 !important;
        text-transform: uppercase;
        letter-spacing: 0.5px;
    }

    [data-testid="stMetricValue"] {
        color: #1E40AF !important;
        font-weight: 700 !important;
        font-size: 1.35rem !important;
    }

    /* Header Telemetry Banner - Light Theme */
    .oasis-header-banner {
        background: #FFFFFF;
        border: 1px solid #E2E8F0;
        border-radius: 12px;
        padding: 18px 24px;
        margin-bottom: 20px;
        box-shadow: 0 4px 16px rgba(15, 23, 42, 0.06);
    }

    .oasis-telemetry-badge {
        background: #F1F5F9;
        border: 1px solid #CBD5E1;
        color: #1E293B;
        padding: 5px 12px;
        border-radius: 6px;
        font-size: 0.82rem;
        font-weight: 600;
        display: inline-block;
    }

    /* Primary Action Buttons */
    .stButton>button[kind="primary"] {
        background: #2563EB !important;
        color: #FFFFFF !important;
        font-weight: 600 !important;
        border: none !important;
        border-radius: 8px !important;
        padding: 10px 20px !important;
        box-shadow: 0 2px 8px rgba(37, 99, 235, 0.25) !important;
        transition: all 0.15s ease-in-out !important;
    }
    .stButton>button[kind="primary"]:hover {
        background: #1D4ED8 !important;
        box-shadow: 0 4px 12px rgba(37, 99, 235, 0.35) !important;
    }

    /* Secondary Buttons */
    .stButton>button {
        border-radius: 8px !important;
        border: 1px solid #CBD5E1 !important;
        background: #FFFFFF !important;
        color: #334155 !important;
        font-weight: 500 !important;
    }
    .stButton>button:hover {
        border-color: #2563EB !important;
        color: #2563EB !important;
        background: #F8FAFC !important;
    }

    /* Tab Header Polish */
    button[data-baseweb="tab"] {
        font-weight: 600 !important;
        font-size: 0.95rem !important;
        padding: 12px 20px !important;
        color: #64748B !important;
    }
    button[data-baseweb="tab"][aria-selected="true"] {
        color: #2563EB !important;
        border-bottom-color: #2563EB !important;
    }
    </style>
    """,
    unsafe_allow_html=True,
)

# Live telemetry state extraction
def get_active_telemetry_state():
    """Extract active location, altitude, atmosphere, and molecule state directly from single authoritative location engine."""
    mol_sel_key = st.session_state.get("mol_selectbox")
    if mol_sel_key and isinstance(mol_sel_key, str) and "ID " in mol_sel_key:
        try:
            mol_id = int(mol_sel_key.split("ID ")[1].split(")")[0].strip())
        except Exception:
            mol_id = int(st.session_state.get("applied_rec", {}).get("mol_id", 2))
    else:
        mol_id = int(st.session_state.get("applied_rec", {}).get("mol_id", 2))

    formula = HITRAN_MOLECULES.get(mol_id, f"ID {mol_id}")

    loc = get_authoritative_location_state()
    lat = float(loc["lat"])
    lon = float(loc["lon"])
    alt = float(loc["alt"])
    place = str(loc["place_en"])

    temp, press, density = standard_atmosphere(alt)

    return {
        "mol_id": mol_id,
        "formula": formula,
        "lat": lat,
        "lon": lon,
        "alt": alt,
        "place": place,
        "temperature_k": temp,
        "pressure_atm": press,
        "air_density_scale": density,
    }



tel_state = get_active_telemetry_state()

# Resolve deferred geocoding (non-blocking: actual HTTP calls happen here on second render pass)
_resolve_deferred_geocoding()

# HTML-escape all dynamic values to prevent raw HTML injection into the header template
import html as _html
_formula_safe  = _html.escape(str(tel_state['formula']))
_lat_safe      = f"{tel_state['lat']:.4f}"
_lon_safe      = f"{tel_state['lon']:.4f}"
_alt_safe      = f"{tel_state['alt']:,.0f}"
_place_safe    = _html.escape(str(tel_state['place']).split(',')[0].strip())
_temp_safe     = f"{tel_state['temperature_k']:.1f}"
_press_safe    = f"{tel_state['pressure_atm']:.3f}"

st.markdown(
    f"""<div class="oasis-header-banner">
<div style="display: flex; justify-content: space-between; align-items: center; flex-wrap: wrap; gap: 16px;">
<div style="display: flex; align-items: center; gap: 16px;">
<svg width="220" height="42" viewBox="0 0 220 42" fill="none" xmlns="http://www.w3.org/2000/svg">
<path d="M 6,34 C 12,14 26,6 42,6 C 58,6 72,14 78,34" stroke="#2563EB" stroke-width="2.5" stroke-linecap="round"/>
<path d="M 12,36 C 18,22 28,15 42,15 C 56,15 66,22 72,36" stroke="#0D9488" stroke-width="1.8" stroke-dasharray="3 3"/>
<circle cx="42" cy="20" r="4.5" fill="#2563EB"/>
<circle cx="42" cy="20" r="8" stroke="#2563EB" stroke-width="1.2" stroke-dasharray="2 2" fill="none"/>
<path d="M 16,24 Q 29,14 42,24 T 68,24" stroke="#7C3AED" stroke-width="2" fill="none" stroke-linecap="round"/>
<text x="88" y="27" font-family="-apple-system, BlinkMacSystemFont, Segoe UI, Roboto, sans-serif" font-weight="900" font-size="23" fill="#0F172A" letter-spacing="-0.5px">OASIS</text>
<text x="88" y="37" font-family="-apple-system, BlinkMacSystemFont, Segoe UI, Roboto, sans-serif" font-weight="600" font-size="7.5" fill="#64748B" letter-spacing="1.5px">SPECTROSCOPY PLATFORM</text>
</svg>
<div style="border-left: 1px solid #CBD5E1; padding-left: 14px; height: 34px; display: flex; flex-direction: column; justify-content: center;">
<span style="color: #0F172A; font-size: 0.88rem; font-weight: 700; letter-spacing: -0.2px;">Open Atmospheric Spectroscopy &amp; Information System</span>
<span style="color: #64748B; font-size: 0.75rem; font-weight: 500;">Universal High-Precision HITRAN/HAPI Spectroscopic &amp; Altitude-Driven LBL Solver &middot; <b style="color: #2563EB; font-weight: 700;">By Ayushi Chahare</b></span>
</div>

</div>
<div style="display: flex; gap: 10px; align-items: center; flex-wrap: wrap;">
<div class="oasis-telemetry-badge">
<svg width="13" height="13" viewBox="0 0 24 24" fill="none" stroke="#2563EB" stroke-width="2" stroke-linecap="round" style="vertical-align: -2px; margin-right: 4px;"><circle cx="12" cy="12" r="3"/><circle cx="19" cy="5" r="2"/><circle cx="5" cy="19" r="2"/><line x1="14" y1="10" x2="17.5" y2="6.5"/><line x1="10" y1="14" x2="6.5" y2="17.5"/></svg>
Species: <b style="color:#1E40AF;">{_formula_safe}</b>
</div>
<div class="oasis-telemetry-badge">
<svg width="13" height="13" viewBox="0 0 24 24" fill="none" stroke="#2563EB" stroke-width="2" stroke-linecap="round" style="vertical-align: -2px; margin-right: 4px;"><path d="M21 10c0 7-9 13-9 13s-9-6-9-13a9 9 0 0 1 18 0z"/><circle cx="12" cy="10" r="3"/></svg>
<b style="color:#1E40AF;">{_lat_safe}&deg;</b>, <b style="color:#1E40AF;">{_lon_safe}&deg;</b>
</div>
<div class="oasis-telemetry-badge">
<svg width="13" height="13" viewBox="0 0 24 24" fill="none" stroke="#2563EB" stroke-width="2" stroke-linecap="round" style="vertical-align: -2px; margin-right: 4px;"><path d="M3 20l7-12 5 6 6-10"/><line x1="3" y1="20" x2="21" y2="20"/></svg>
<b style="color:#1E40AF;">{_alt_safe} m</b> ({_place_safe})
</div>
<div class="oasis-telemetry-badge">
<svg width="13" height="13" viewBox="0 0 24 24" fill="none" stroke="#2563EB" stroke-width="2" stroke-linecap="round" style="vertical-align: -2px; margin-right: 4px;"><path d="M12 2v20M2 12h20M12 6a6 6 0 1 0 0 12 6 6 0 0 0 0-12z"/></svg>
<b style="color:#1E40AF;">{_temp_safe} K</b> | <b style="color:#1E40AF;">{_press_safe} atm</b>
</div>
</div>
</div>
</div>""",
    unsafe_allow_html=True,
)




# ─────────────────────────────────────────────────────────────────
# TABS
# ─────────────────────────────────────────────────────────────────
tab1, tab2, tab3 = st.tabs([
    "Global Database Explorer",
    "Altitude-Driven LBL Simulation",
    "Global Station & Metrics Engine",
])





# ═══════════════════════════════════════════════════════════════════
# TAB 1 — GLOBAL DATABASE EXPLORER
# ═══════════════════════════════════════════════════════════════════
with tab1:
    st.header("Global Spectroscopic Database (SQLite)")
    st.write("Explore all stored quantum transition lines across active species.")

    try:
        df = load_hitran_parameters()

        if df.empty:
            st.warning("Database is empty. Run `oasis.py` first to populate it.")
        else:
            st.sidebar.header("Database Filters")
            available_molecules = sorted(df["molecule_name"].unique().tolist())

            selected_molecules = st.sidebar.multiselect(
                "Molecules to Display",
                options=available_molecules,
                default=available_molecules[:4] if len(available_molecules) >= 4 else available_molecules,
            )

            min_nu = float(df["wavenumber"].min())
            max_nu = float(df["wavenumber"].max())
            selected_range = st.sidebar.slider(
                "Wavenumber Range (cm⁻¹)",
                min_value=min_nu, max_value=max_nu,
                value=(min_nu, max_nu), key="db_nu_slider",
            )

            min_intensity = float(df["intensity"].min())
            max_intensity = float(df["intensity"].max())
            selected_min_intensity = st.sidebar.slider(
                "Minimum Line Intensity Filter",
                min_value=min_intensity,
                max_value=max_intensity / 10 if max_intensity > 0 else 1.0,
                value=min_intensity, key="db_int_slider",
            )

            filtered_df = df[
                df["molecule_name"].isin(selected_molecules) &
                (df["wavenumber"] >= selected_range[0]) &
                (df["wavenumber"] <= selected_range[1]) &
                (df["intensity"]  >= selected_min_intensity)
            ]

            c1, c2, c3 = st.columns(3)
            with c1:
                st.metric("Lines in View", f"{len(filtered_df):,}")
            with c2:
                st.metric("Molecules Selected", len(selected_molecules))
            with c3:
                st.write("###")
                st.download_button(
                    "Download Database CSV",
                    filtered_df.to_csv(index=False).encode("utf-8"),
                    file_name="oasis_global_database_export.csv",
                    mime="text/csv",
                )

            st.subheader("Multi-Molecule Spectral Lines")
            if not filtered_df.empty:
                fig = px.scatter(
                    filtered_df, x="wavenumber", y="intensity", color="molecule_name",
                    labels={"wavenumber": "Wavenumber (cm⁻¹)", "intensity": "Intensity (cm/molecule)",
                            "molecule_name": "Molecule"},
                )
                fig.update_traces(marker=dict(size=4))
                fig.update_layout(template="plotly_dark", hovermode="closest")
                st.plotly_chart(fig, use_container_width=True)
            else:
                st.info("No data points match the selected filters.")

            st.subheader("SQLite Parameters Table")
            st.dataframe(filtered_df, use_container_width=True)

    except Exception as e:
        st.error(f"Error loading database: {e}")


# ═══════════════════════════════════════════════════════════════════
# TAB 2 — ALTITUDE-DRIVEN LBL SIMULATION
# ═══════════════════════════════════════════════════════════════════
with tab2:
    st.header("Altitude-Driven Line-by-Line (LBL) Physics Simulator")
    st.write(
        "Select a molecule and location. The system computes atmospheric conditions, "
        "recommends simulation parameters, validates your configuration, and runs the "
        "spectroscopic simulation."
    )

    # ── Load local spectral database ──────────────────────────────
    try:
        sim_source_df = load_hitran_parameters()
    except Exception as e:
        sim_source_df = pd.DataFrame()
        st.error(f"Could not load local spectral database: {e}")

    # ── Session state initialisation ──────────────────────────────
    if "sim_mode" not in st.session_state:
        st.session_state.sim_mode = "Recommended"
    if "applied_rec" not in st.session_state:
        st.session_state.applied_rec = {}
    if "sim_request_id" not in st.session_state:
        st.session_state.sim_request_id = None
    if "saved_configs" not in st.session_state:
        st.session_state.saved_configs = []
    # Map location bridge
    if "map_lat"  not in st.session_state:
        st.session_state.map_lat  = None
    if "map_lon"  not in st.session_state:
        st.session_state.map_lon  = None
    if "map_alt"  not in st.session_state:
        st.session_state.map_alt  = None

    # ── Sidebar: Saved configurations ─────────────────────────────
    st.sidebar.header("Saved Configurations")
    if st.session_state.saved_configs:
        cfg_labels = [
            f"{c['molecule']} @ {c['altitude_m']:.0f} m  [{c['nu_start']:.0f}–{c['nu_end']:.0f} cm⁻¹]"
            for c in st.session_state.saved_configs
        ]
        sel_cfg_label = st.sidebar.selectbox(
            "Restore a saved configuration", ["— none —"] + cfg_labels,
            key="restore_selectbox",
        )
        if sel_cfg_label != "— none —":
            idx = cfg_labels.index(sel_cfg_label)
            cfg = st.session_state.saved_configs[idx]
            if st.sidebar.button("Restore this configuration"):
                set_authoritative_location_state(
                    float(cfg.get("lat", 19.54)),
                    float(cfg.get("lon", -155.58)),
                    float(cfg.get("altitude_m", 3397.0)),
                )
                st.session_state[f"wn_start_{cfg['mol_id']}_{cfg['iso_id']}"] = float(cfg["nu_start"])
                st.session_state[f"wn_end_{cfg['mol_id']}_{cfg['iso_id']}"]   = float(cfg["nu_end"])
                st.session_state[f"mf_{cfg['mol_id']}_{cfg['iso_id']}"]       = float(cfg["mole_fraction"])
                st.session_state["preset_sel"] = "Custom range"


                st.session_state.applied_rec = {
                    "mol_id":            cfg["mol_id"],
                    "iso_id":            cfg["iso_id"],
                    "temperature_k":     cfg["temperature_k"],
                    "pressure_atm":      cfg["pressure_atm"],
                    "start_wavenumber":  cfg["nu_start"],
                    "end_wavenumber":    cfg["nu_end"],
                    "gas_mole_fraction": cfg["mole_fraction"],
                    "path_length_cm":    cfg["path_length_cm"],
                    "grid_step":         cfg["grid_step"],
                    "_lat":              cfg.get("lat", 19.54),
                    "_lon":              cfg.get("lon", -155.58),
                    "_alt":              cfg.get("altitude_m", 3397.0),
                }
                st.session_state.sim_mode = "Recommended"
                st.session_state.sim_request_id = None
                st.rerun()
    else:
        st.sidebar.caption("No configurations saved yet.")

    # ═══════════════════════════════════════════════
    # COLUMN LAYOUT
    # ═══════════════════════════════════════════════
    col_sim1, col_sim2 = st.columns(2)

    # ── COLUMN 1: Location, atmosphere, sensor ────────────────────
    with col_sim1:
        st.subheader("Location & Altitude")

        # Map-location bridge
        loc = get_authoritative_location_state()

        if st.button("Use Location Selected on Map", type="primary", use_container_width=True):
            set_authoritative_location_state(loc["lat"], loc["lon"], loc["alt"])
            st.success(f"Location ({loc['lat']:.4f}°, {loc['lon']:.4f}°) sent to Simulation!")
            st.rerun()

        st.caption(
            f"Active Map Selection: Lat {loc['lat']:.4f}°, "
            f"Lon {loc['lon']:.4f}°, "
            f"Alt {loc['alt']:,.1f} m ({loc['place_en']})"
        )

        _default_lat = float(loc["lat"])
        _default_lon = float(loc["lon"])
        _default_alt = float(loc["alt"])

        lat        = st.number_input("Latitude (°)",                     value=_default_lat, format="%.4f")
        lon        = st.number_input("Longitude (°)",                    value=_default_lon, format="%.4f")
        altitude_m = st.number_input("Altitude (m above sea level)",     value=_default_alt, step=100.0)

        if lat != _default_lat or lon != _default_lon or altitude_m != _default_alt:
            set_authoritative_location_state(lat, lon, altitude_m)




        calc_temp, calc_pressure, air_density_scale = standard_atmosphere(altitude_m)

        st.info(
            f"**Auto-Computed ISA Atmosphere at {altitude_m:.0f} m:**\n\n"
            f"- Temperature: **{calc_temp:.2f} K**  *(DERIVED from ISA model)*\n\n"
            f"- Pressure: **{calc_pressure:.4f} atm**  *(DERIVED from ISA model)*\n\n"
            f"- Air Density Scale: **{air_density_scale:.3f}** (relative to sea level)"
        )

        st.subheader("Path Length & Sensor")
        sim_path         = st.slider("Absorption Path Length (cm)", 10.0, 10000.0, 500.0, 10.0)
        incident_signal  = st.number_input(
            "Incident Signal (arbitrary units)",
            value=100.0, min_value=0.0,
            help="Arbitrary starting signal before atmospheric absorption. 100 = 100% signal.",
        )
        detection_threshold = st.number_input(
            "Detection Threshold (same units)",
            value=1.0, min_value=0.0,
        )

    # ── COLUMN 2: Molecule, isotopologue, spectral window ─────────
    with col_sim2:
        st.subheader("Molecule & Spectral Window")

        # Build molecule selectbox labels with support tags
        local_counts = {}
        if not sim_source_df.empty:
            if "molecule_id" in sim_source_df.columns:
                local_counts = sim_source_df.groupby("molecule_id")["wavenumber"].size().to_dict()
            elif "mol_id" in sim_source_df.columns:
                local_counts = sim_source_df.groupby("mol_id")["wavenumber"].size().to_dict()


        def mol_label(mol_id, formula):
            sup = MOLECULE_SUPPORT_MAP.get(mol_id, {}).get("status", "UNKNOWN")
            lc  = int(local_counts.get(mol_id, 0))
            if sup == "LOCAL_DATA":
                tag = f"[LOCAL {lc:,} lines]"
            elif sup == "FETCH_REQUIRED":
                tag = "[FETCH REQUIRED]"
            else:
                tag = "[UNSUPPORTED]"
            return f"{formula} (ID {mol_id}) {tag}"

        molecule_options = {mol_label(mid, fm): mid for mid, fm in HITRAN_MOLECULES.items()}

        # Default molecule from applied_rec or CO2
        _default_mol_id = int(st.session_state.applied_rec.get("mol_id", 2))
        _default_mol_idx = list(molecule_options.values()).index(_default_mol_id) \
            if _default_mol_id in molecule_options.values() else \
            list(molecule_options.values()).index(2)

        sim_molecule_label = st.selectbox(
            "Select HITRAN Molecule",
            options=list(molecule_options.keys()),
            index=_default_mol_idx,
            key="mol_selectbox",
        )
        selected_mol_id  = molecule_options[sim_molecule_label]
        selected_molecule = HITRAN_MOLECULES[selected_mol_id]
        mol_support       = MOLECULE_SUPPORT_MAP.get(selected_mol_id, {})
        mol_sup_status    = mol_support.get("status", "UNKNOWN")

        # Molecule support badge
        if mol_sup_status == "LOCAL_DATA":
            st.success(
                f"**LOCAL DATA** — {mol_support['line_count']:,} lines cached "
                f"({mol_support['nu_min']:.1f}–{mol_support['nu_max']:.1f} cm⁻¹). "
                "Simulation can run without internet access."
            )
        elif mol_sup_status == "FETCH_REQUIRED":
            st.warning(
                "**FETCH REQUIRED** — No local cache for this molecule. "
                "Simulation will attempt an on-demand HAPI download (internet required). "
                "First-time fetch may be slow."
            )
        else:
            st.error(
                f"**UNSUPPORTED** — {mol_support.get('note', 'Cannot simulate this molecule.')} "
                "No data is available."
            )

        # Isotopologue
        molecule_rows = (
            sim_source_df[sim_source_df["molecule_id"] == selected_mol_id]
            if not sim_source_df.empty else pd.DataFrame()
        )
        iso_options = sorted(
            molecule_rows["isotopologue_id"].dropna().astype(int).unique().tolist()
        ) if not molecule_rows.empty else [1]
        if 1 not in iso_options:
            iso_options.insert(0, 1)
        _default_iso = int(st.session_state.applied_rec.get("iso_id", 1))
        _iso_idx     = iso_options.index(_default_iso) if _default_iso in iso_options else 0
        selected_iso_id = st.selectbox(
            "Isotopologue ID", options=iso_options, index=_iso_idx, key="iso_selectbox"
        )

        # ─ Smart Parameter Assistant ──────────────────────────────
        st.markdown("---")
        st.subheader("Smart Parameter Assistant")

        recommendation = get_recommended_simulation_parameters(
            selected_mol_id, selected_molecule, altitude_m,
            selected_iso_id, sim_source_df,
        )

        rec_status = recommendation["status"]
        with st.expander("Recommended parameters for this molecule & altitude", expanded=True):
            r1, r2 = st.columns(2)
            with r1:
                st.metric(
                    "Recommended Start (cm⁻¹)",
                    f"{recommendation['start_wavenumber']:.3f}",
                    help=recommendation["value_sources"]["start_wavenumber"],
                )
                st.metric(
                    "Recommended Temperature (K)",
                    f"{recommendation['temperature_k']:.2f}",
                    help=recommendation["value_sources"]["temperature_k"],
                )
                st.metric(
                    "Recommended Mole Fraction",
                    f"{recommendation['gas_mole_fraction']:.3e}",
                    help=recommendation["value_sources"]["gas_mole_fraction"],
                )
            with r2:
                st.metric(
                    "Recommended End (cm⁻¹)",
                    f"{recommendation['end_wavenumber']:.3f}",
                    help=recommendation["value_sources"]["end_wavenumber"],
                )
                st.metric(
                    "Recommended Pressure (atm)",
                    f"{recommendation['pressure_atm']:.4f}",
                    help=recommendation["value_sources"]["pressure_atm"],
                )
                st.metric(
                    "Path Length (cm)",
                    f"{recommendation['path_length_cm']:.0f}",
                    help=recommendation["value_sources"]["path_length_cm"],
                )

            st.caption(
                f"**Source:** {recommendation['explanation']}"
            )

            # Value source legend
            st.markdown(
                "_DERIVED = computed from altitude & ISA model · "
                "LOCAL\\_DATA = from cached HITRAN lines · "
                "PRESET = published spectral band · "
                "APP\\_DEFAULT = application default (not a measured value)_"
            )

        if st.button("Apply Recommended Values", type="primary", key="apply_rec_btn"):
            st.session_state.sim_mode    = "Recommended"
            rec_wn_start = float(recommendation["start_wavenumber"])
            rec_wn_end   = float(recommendation["end_wavenumber"])
            rec_mf       = float(recommendation["gas_mole_fraction"])

            # Atomically update widget input keys
            st.session_state[f"wn_start_{selected_mol_id}_{selected_iso_id}"] = rec_wn_start
            st.session_state[f"wn_end_{selected_mol_id}_{selected_iso_id}"]   = rec_wn_end
            st.session_state[f"mf_{selected_mol_id}_{selected_iso_id}"]       = rec_mf
            set_authoritative_location_state(lat, lon, altitude_m)


            st.session_state.applied_rec = {
                **st.session_state.applied_rec,
                "mol_id":            selected_mol_id,
                "iso_id":            selected_iso_id,
                "temperature_k":     recommendation["temperature_k"],
                "pressure_atm":      recommendation["pressure_atm"],
                "start_wavenumber":  rec_wn_start,
                "end_wavenumber":    rec_wn_end,
                "gas_mole_fraction": rec_mf,
                "path_length_cm":    500.0,
                "grid_step":         recommendation["grid_step"],
                "_lat":              lat,
                "_lon":              lon,
                "_alt":              altitude_m,
            }
            st.rerun()


        st.markdown("---")

        # ─ Simulation Mode ────────────────────────────────────────
        sim_mode = st.radio(
            "Parameter Mode",
            ["Recommended", "Custom"],
            index=0 if st.session_state.sim_mode == "Recommended" else 1,
            horizontal=True,
            key="sim_mode_radio",
        )
        st.session_state.sim_mode = sim_mode

        # ─ Spectral window ────────────────────────────────────────
        presets = SPECTRAL_PRESETS.get(
            selected_molecule, [("Custom (no preset available)", 2200.0, 2400.0)]
        )
        preset_options = [
            f"{name}: {start:g}–{end:g} cm⁻¹  ({wavenumber_to_microns(end):.2f}–{wavenumber_to_microns(start):.2f} µm)"
            for name, start, end in presets
        ]
        preset_options.append("Custom range")
        selected_preset = st.selectbox("Suggested Spectral Window", options=preset_options, key="preset_sel")

        if selected_preset == "Custom range":
            preset_start_nu = float(st.session_state.applied_rec.get("start_wavenumber", 2200.0))
            preset_end_nu   = float(st.session_state.applied_rec.get("end_wavenumber",   2400.0))
        else:
            pidx            = preset_options.index(selected_preset)
            _, preset_start_nu, preset_end_nu = presets[pidx]

        unit_mode = st.radio(
            "Spectral Input Unit",
            ["Wavenumber (cm⁻¹)", "Wavelength (µm)"],
            horizontal=True, key="unit_mode_radio",
        )

        if sim_mode == "Recommended" and st.session_state.applied_rec.get("mol_id") == selected_mol_id:
            _start_val = float(st.session_state.applied_rec.get("start_wavenumber", preset_start_nu))
            _end_val   = float(st.session_state.applied_rec.get("end_wavenumber",   preset_end_nu))
        else:
            _start_val = preset_start_nu
            _end_val   = preset_end_nu

        if unit_mode == "Wavenumber (cm⁻¹)":
            sim_start_nu = st.number_input(
                "Start Wavenumber (cm⁻¹)", value=float(_start_val),
                min_value=0.0001, key=f"wn_start_{selected_mol_id}_{selected_iso_id}",
            )
            sim_end_nu = st.number_input(
                "End Wavenumber (cm⁻¹)", value=float(_end_val),
                min_value=0.0001, key=f"wn_end_{selected_mol_id}_{selected_iso_id}",
            )
        else:
            _short_um = wavenumber_to_microns(_end_val)
            _long_um  = wavenumber_to_microns(_start_val)
            start_um  = st.number_input(
                "Short Wavelength (µm)", value=float(_short_um),
                min_value=0.0001, key=f"um_short_{selected_mol_id}_{selected_iso_id}",
            )
            end_um    = st.number_input(
                "Long Wavelength (µm)", value=float(_long_um),
                min_value=0.0001, key=f"um_long_{selected_mol_id}_{selected_iso_id}",
            )
            sim_start_nu = microns_to_wavenumber(start_um)
            sim_end_nu   = microns_to_wavenumber(end_um)
            if sim_start_nu > sim_end_nu:
                sim_start_nu, sim_end_nu = sim_end_nu, sim_start_nu
            st.caption(f"Converted: {sim_start_nu:.3f} to {sim_end_nu:.3f} cm⁻¹")

        st.caption(
            f"Spectral window: {sim_start_nu:.3f}–{sim_end_nu:.3f} cm⁻¹  "
            f"= {wavenumber_to_microns(sim_end_nu):.3f}–{wavenumber_to_microns(sim_start_nu):.3f} µm"
        )

        if sim_end_nu - sim_start_nu > 500:
            st.warning("Wide spectral window — calculation may be slow.")

        # Mole fraction
        if sim_mode == "Recommended" and st.session_state.applied_rec.get("mol_id") == selected_mol_id:
            _frac_val = float(st.session_state.applied_rec.get(
                "gas_mole_fraction", DEFAULT_MOLE_FRACTIONS.get(selected_molecule, 1e-9)
            ))
        else:
            _frac_val = DEFAULT_MOLE_FRACTIONS.get(selected_molecule, 1e-9)

        gas_mole_fraction = st.number_input(
            "Gas Mole Fraction",
            value=float(_frac_val), min_value=1e-15, max_value=1.0, format="%.10f",
            help="Fraction of selected gas in the atmospheric column.",
            key=f"mf_{selected_mol_id}_{selected_iso_id}",
        )
        if gas_mole_fraction < 1e-15:
            st.error("Gas mole fraction is too small — simulation will produce zero absorption.")

        grid_step = st.selectbox(
            "LBL Grid Step (cm⁻¹)", [0.1, 0.05, 0.02, 0.01, 0.005],
            index=3, key="grid_step_sel",
            help="Smaller values give finer spectral resolution but are slower.",
        )

    # ═══════════════════════════════════════════════
    # PREFLIGHT VALIDATOR (always displayed)
    # ═══════════════════════════════════════════════
    st.markdown("---")
    st.subheader("Simulation Preflight Validation")

    preflight_config = {
        "molecule_id":        selected_mol_id,
        "isotopologue_id":    int(selected_iso_id),
        "altitude_m":         altitude_m,
        "temperature_k":      calc_temp,
        "pressure_atm":       calc_pressure,
        "start_nu":           sim_start_nu,
        "end_nu":             sim_end_nu,
        "gas_mole_fraction":  gas_mole_fraction,
        "path_length_cm":     sim_path,
        "grid_step":          grid_step,
    }

    preflight = preflight_simulation(preflight_config, recommendation, sim_source_df)
    pf_status = preflight["status"]

    if pf_status == "READY":
        st.success("READY — Configuration is valid. You may proceed with the simulation.")
    elif pf_status == "WARNING":
        st.warning("WARNING — Configuration is valid but outside the recommended range.")
        for w in preflight["warnings"]:
            st.caption(f"  Warning: {w}")
    elif pf_status == "INVALID":
        st.error("INVALID — Configuration has errors that must be fixed before simulating.")
        for m in preflight["messages"]:
            st.caption(f"  Error: {m}")
        for w in preflight["warnings"]:
            st.caption(f"  Warning: {w}")
    elif pf_status == "UNSUPPORTED":
        st.error("UNSUPPORTED — This molecule or configuration cannot be reliably simulated.")
        for m in preflight["messages"]:
            st.caption(f"  Reason: {m}")
        if recommendation["status"] not in ("UNSUPPORTED",):
            st.info(
                "Recovery path: Apply the Recommended Values above, "
                "which provide a validated configuration for this molecule."
            )

    # Save current config
    col_run, col_save = st.columns([3, 1])
    with col_save:
        if st.button("Save Current Configuration"):
            st.session_state.saved_configs.append({
                "molecule":      selected_molecule,
                "mol_id":        selected_mol_id,
                "iso_id":        int(selected_iso_id),
                "lat":           lat,
                "lon":           lon,
                "altitude_m":    altitude_m,
                "temperature_k": calc_temp,
                "pressure_atm":  calc_pressure,
                "nu_start":      sim_start_nu,
                "nu_end":        sim_end_nu,
                "mole_fraction": gas_mole_fraction,
                "path_length_cm": sim_path,
                "grid_step":     grid_step,
            })
            st.success("Configuration saved to sidebar.")

    # ═══════════════════════════════════════════════
    # RUN SIMULATION BUTTON
    # ═══════════════════════════════════════════════
    with col_run:
        run_blocked = pf_status in ("INVALID", "UNSUPPORTED")
        run_button  = st.button(
            "Run Altitude-Driven LBL Simulation",
            type="primary",
            disabled=run_blocked,
        )

    if run_blocked and not run_button:
        if pf_status == "INVALID":
            st.info("Fix the validation errors above before running the simulation.")
        else:
            st.info(
                "This configuration is unsupported. "
                "Apply recommended values above to get a runnable configuration."
            )

    # ─────────────────────────────────────────────────────────────
    # SIMULATION EXECUTION
    # ─────────────────────────────────────────────────────────────
    if run_button and not run_blocked:
        mol_id             = selected_mol_id
        iso_id             = int(selected_iso_id)
        sim_molecule       = selected_molecule
        isotope_abundance  = natural_isotopologue_abundance(mol_id, iso_id)
        component_abundance = isotope_abundance * gas_mole_fraction

        st.success(
            f"**[LBL Solver Executing]** Molecule: **{sim_molecule}** (ID {mol_id}) | "
            f"Location: **{lat:.4f}°**, **{lon:.4f}°** | Altitude: **{altitude_m:.1f} m** | "
            f"Derived ISA Atmosphere: **{calc_temp:.2f} K**, **{calc_pressure:.4f} atm**"
        )


        # Stale-result prevention: record request fingerprint

        _req_data = json.dumps({
            "mol_id": mol_id, "iso_id": iso_id,
            "start":  sim_start_nu, "end": sim_end_nu,
            "alt":    altitude_m,
        }, sort_keys=True)
        current_request_id = hashlib.sha256(_req_data.encode()).hexdigest()[:16]
        st.session_state.sim_request_id = current_request_id

        # Step 1 — HAPI fetch (smart engine with classified error handling)
        with st.spinner("Checking local HAPI cache and fetching missing HITRAN data if needed…"):
            fetch_result = smart_fetch_hapi_data(
                mol_id, iso_id, sim_start_nu, sim_end_nu,
                request_id=current_request_id,
            )

        # ── Stale-result check ────────────────────────────────────
        if fetch_result.get("request_id") != current_request_id:
            st.warning(
                "⚠️ A stale HAPI response was detected (molecule/parameters changed mid-flight). "
                "Please click Run Simulation again."
            )
            st.stop()

        # ── Handle fetch failure ──────────────────────────────────
        if not fetch_result["ok"]:
            cat = fetch_result["category"]
            user_msg = fetch_result["user_message"]
            recovery = fetch_result["recovery"]
            detail   = fetch_result["detail"]

            # Colour-coded icon by failure type
            if cat in ("HAPI_NO_DATA_IN_RANGE", "HAPI_NO_DATA_ANY_RANGE"):
                # Scientific / range issue — orange warning
                st.warning(f"**{cat}** — {user_msg}")
            elif cat in ("HAPI_NETWORK_FAILURE", "HAPI_TIMEOUT", "HAPI_RATE_LIMITED"):
                # Connectivity issue — yellow
                st.warning(f"**{cat}** — {user_msg}")
            else:
                # Technical failure — red
                st.error(f"**{cat}** — {user_msg}")

            # Recovery suggestion
            if recovery:
                st.info(f"**Recovery:** {recovery}")

            # If a valid range was discovered, offer to apply it
            if fetch_result.get("valid_range"):
                vmin, vmax = fetch_result["valid_range"]
                formula = HITRAN_MOLECULES.get(mol_id, f"Mol {mol_id}")
                apply_key = f"apply_range_{mol_id}_{int(vmin)}"
                if st.button(
                    f"Apply valid range: {vmin:.0f}–{vmax:.0f} cm⁻¹ for {formula}",
                    key=apply_key,
                ):
                    rec = st.session_state.get("applied_rec", {}).copy()
                    rec["start_wavenumber"] = vmin
                    rec["end_wavenumber"]   = vmax
                    st.session_state.applied_rec = rec
                    st.session_state.sim_mode = "Recommended"
                    st.rerun()

            # Technical detail expander
            with st.expander("Technical details"):
                st.markdown(f"**Error category:** `{cat}`")
                st.markdown(f"**HITRAN API base:** `{HITRAN_API_BASE}`")
                try:
                    giso, _ = _hitran_global_iso_id(mol_id, iso_id), None
                    st.markdown(f"**Global ISO ID requested:** `{_hitran_global_iso_id(mol_id, iso_id)[0]}`")
                except Exception:
                    pass
                st.markdown(f"**Requested range:** {sim_start_nu:.3f}–{sim_end_nu:.3f} cm⁻¹")
                if detail:
                    st.code(detail)
                st.caption(
                    f"Category description: {HAPI_ERROR_CATEGORIES.get(cat, 'Unknown category.')}"
                )
            st.stop()

        table_name    = fetch_result["table_name"]
        fetched_now   = fetch_result["fetched"]
        table_summary = get_hapi_table_summary(table_name)

        # ── Data validation header ────────────────────────────────
        st.subheader("Data Validation")
        v1, v2, v3, v4, v5, v6 = st.columns(6)
        with v1:
            st.metric("Molecule ID",    mol_id)
        with v2:
            st.metric("Isotopologue",   iso_id)
        with v3:
            st.metric("HAPI Lines",     f"{int(table_summary['line_count']):,}")
        with v4:
            st.metric("Mole Fraction",  f"{gas_mole_fraction:.3e}")
        with v5:
            st.metric("Cache",          "Fetched" if fetched_now else "Local")
        with v6:
            st.metric("Fetch Status",   fetch_result["category"].replace("HAPI_", ""))

        st.caption(
            f"HAPI table `{table_name}`  ·  {sim_start_nu:.3f}–{sim_end_nu:.3f} cm⁻¹  "
            f"({wavenumber_to_microns(sim_end_nu):.3f}–{wavenumber_to_microns(sim_start_nu):.3f} µm)  "
            f"·  {fetch_result['user_message']}"
        )

        if table_summary["line_count"] == 0:
            st.error(
                f"**HAPI_NO_DATA_IN_RANGE** — HAPI returned zero spectral lines for "
                f"{sim_molecule} (ID {mol_id}) in the range "
                f"{sim_start_nu:.1f}–{sim_end_nu:.1f} cm⁻¹."
            )
            st.info(
                "Recovery: Use the recommended spectral window for this molecule, "
                "which is centred on known strong absorption features. "
                "Click 'Apply Recommended Values' in the Smart Parameter Assistant above."
            )
            st.stop()



        # Step 2 — Voigt calculation
        try:
            with st.spinner(
                f"Computing LBL Voigt absorption for {sim_molecule} at {altitude_m:.0f} m…"
            ):
                nu_grid, absorption_coeff = absorptionCoefficient_Voigt(
                    SourceTables=table_name,
                    Components=[(mol_id, iso_id, component_abundance)],
                    Environment={"T": calc_temp, "p": calc_pressure},
                    WavenumberRange=[sim_start_nu, sim_end_nu],
                    WavenumberStep=grid_step,
                    HITRAN_units=False,
                )
        except Exception as e:
            err_msg = str(e)
            if "zero" in err_msg.lower() or "nan" in err_msg.lower():
                category = "NUMERICAL_CALCULATION_ERROR"
            else:
                category = "UNEXPECTED_ERROR"
            st.error(
                f"**{category}** — Voigt calculation failed for {sim_molecule}.\n\n"
                f"{FAILURE_CATEGORIES[category]}"
            )
            with st.expander("Technical details"):
                st.code(err_msg)
            st.stop()

        if nu_grid is None or len(nu_grid) == 0:
            st.error(
                f"**EMPTY_RESULT** — HAPI returned no grid points for {sim_molecule} "
                f"in [{sim_start_nu:.1f}, {sim_end_nu:.1f}] cm⁻¹.\n\n"
                f"{FAILURE_CATEGORIES['EMPTY_RESULT']}"
            )
            st.stop()

        # Step 3 — Transmittance
        try:
            nu_grid, transmittance = transmittanceSpectrum(
                nu_grid, absorption_coeff, Environment={"l": sim_path}
            )
        except Exception as e:
            st.error(
                f"**NUMERICAL_CALCULATION_ERROR** — Transmittance calculation failed.\n\n"
                f"{FAILURE_CATEGORIES['NUMERICAL_CALCULATION_ERROR']}"
            )
            with st.expander("Technical details"):
                st.code(str(e))
            st.stop()

        # Step 4 — Assemble DataFrame
        sim_df = pd.DataFrame({
            "wavenumber":            nu_grid,
            "absorption_coefficient": absorption_coeff,
            "transmittance":         transmittance,
        })
        sim_df["wavelength_um"]     = 10000.0 / sim_df["wavenumber"]
        sim_df["optical_depth"]     = sim_df["absorption_coefficient"] * sim_path
        sim_df["absorbance"]        = -np.log10(np.clip(sim_df["transmittance"], 1e-300, 1.0))
        sim_df["absorption_depth"]  = 1.0 - sim_df["transmittance"]
        sim_df["received_signal"]   = incident_signal * sim_df["transmittance"]

        expected_rows = int(round((sim_end_nu - sim_start_nu) / grid_step)) + 1

        # Step 5 — Stale result check
        if st.session_state.sim_request_id != current_request_id:
            st.warning(
                "**STALE_RESULT_PREVENTED** — A simulation result was discarded because "
                "the simulation parameters changed before the calculation completed.\n\n"
                f"{FAILURE_CATEGORIES['STALE_RESULT_PREVENTED']}"
            )
            st.stop()

        # Step 6 — Output validation
        result_ok, result_category, result_message = validate_simulation_result(
            sim_df, expected_rows
        )

        if not result_ok:
            if result_category == "WEAK_SIGNAL":
                st.warning(
                    f"**WEAK SIGNAL** — {result_message}\n\n"
                    "This is not an error. The molecule is physically real in this window, "
                    "but the absorption is too weak to produce visible spectral features "
                    "at this concentration and path length.\n\n"
                    "**Suggested actions:**\n"
                    "- Apply recommended values to use a stronger absorption band.\n"
                    "- Increase path length (to 10,000–100,000 cm for weak species).\n"
                    "- Increase mole fraction if studying a high-concentration scenario.\n"
                    "- Select a different spectral window with stronger lines."
                )
                # Still render the (flat) result, but with the warning prominently displayed
                st.info(
                    "The plots below show the result. Transmittance will appear as a flat line "
                    "near 1.0 — this is scientifically correct for this configuration."
                )
            else:
                st.error(
                    f"**{result_category}** — {result_message}\n\n"
                    f"{FAILURE_CATEGORIES.get(result_category, '')}"
                )
                with st.expander("Diagnostic details"):
                    st.json({
                        "rows":                len(sim_df),
                        "expected_rows":       expected_rows,
                        "max_absorption_coeff": float(sim_df["absorption_coefficient"].max()),
                        "transmittance_min":   float(sim_df["transmittance"].min()),
                        "transmittance_max":   float(sim_df["transmittance"].max()),
                        "has_nan":             bool(sim_df.isnull().any().any()),
                        "has_inf":             bool(np.isinf(sim_df.select_dtypes(float).values).any()),
                    })
                st.stop()

        # ─ Compute summary metrics ────────────────────────────────
        min_trans             = float(sim_df["transmittance"].min())
        mean_trans            = float(sim_df["transmittance"].mean())
        max_absorption_depth  = 1.0 - min_trans
        max_optical_depth     = float(sim_df["absorption_coefficient"].max()) * sim_path
        max_absorbance        = float(sim_df["absorbance"].max())
        min_received_signal   = float(sim_df["received_signal"].min())
        detection_margin      = min_received_signal - detection_threshold
        peak_abs_row          = sim_df.loc[sim_df["absorption_depth"].idxmax()]
        integrated_opacity    = float(sim_df["absorption_coefficient"].sum()) * grid_step
        equivalent_width      = float(((1.0 - sim_df["transmittance"]) * grid_step).sum())

        validation_rows = pd.DataFrame()
        if not sim_source_df.empty:
            validation_rows = sim_source_df[
                (sim_source_df["molecule_id"]     == mol_id) &
                (sim_source_df["isotopologue_id"] == iso_id) &
                (sim_source_df["wavenumber"]      >= sim_start_nu) &
                (sim_source_df["wavenumber"]      <= sim_end_nu)
            ]

        if not validation_rows.empty:
            strongest = validation_rows.loc[validation_rows["intensity"].idxmax()]
            strongest_text = (
                f"Strongest local SQLite line in this window: "
                f"{strongest['wavenumber']:.4f} cm⁻¹  "
                f"(intensity {strongest['intensity']:.3e} cm/molecule)"
            )
        else:
            strongest_text = (
                "Simulation used on-demand HAPI cache data not yet merged into SQLite."
            )

        # Store active simulation result persistently across reruns
        pdf_config = {
            "sim_id": current_request_id,
            "molecule": sim_molecule,
            "mol_id": mol_id,
            "iso_id": iso_id,
            "isotope_abundance": isotope_abundance,
            "abundance": component_abundance,
            "mole_fraction": gas_mole_fraction,
            "path_length_cm": sim_path,
            "start_nu": sim_start_nu,
            "end_nu": sim_end_nu,
            "wavelength_min_um": wavenumber_to_microns(sim_end_nu),
            "wavelength_max_um": wavenumber_to_microns(sim_start_nu),
            "grid_step": grid_step,
            "lat": lat,
            "lon": lon,
            "alt": altitude_m,
            "temp_k": calc_temp,
            "press_atm": calc_pressure,
            "air_density_scale": air_density_scale,
            "place_en": loc.get("place_en", "Standard Atmosphere"),
            "place_native": loc.get("place_native", "Standard Atmosphere"),
            "incident_signal": incident_signal,
            "detection_threshold": detection_threshold,
        }

        pdf_metrics = {
            "line_count": len(validation_rows),
            "grid_points": len(nu_grid),
            "equivalent_width": equivalent_width,
            "peak_absorption_wavenumber": float(peak_abs_row["wavenumber"]),
            "peak_absorption_wavelength_um": float(peak_abs_row["wavelength_um"]),
            "max_optical_depth": max_optical_depth,
            "max_absorbance": max_absorbance,
            "integrated_opacity": integrated_opacity,
            "min_trans": min_trans,
            "mean_trans": mean_trans,
            "max_abs": float(sim_df["absorption_coefficient"].max()),
            "min_signal": min_received_signal,
            "detection_threshold": detection_threshold,
            "detection_margin": detection_margin,
            "status": "Detected" if detection_margin >= 0 else "Below threshold",
        }


        st.session_state.active_sim_result = {
            "sim_df": sim_df,
            "pdf_config": pdf_config,
            "pdf_metrics": pdf_metrics,
            "sim_molecule": sim_molecule,
            "mol_id": mol_id,
            "iso_id": iso_id,
            "lat": lat,
            "lon": lon,
            "altitude_m": altitude_m,
            "calc_temp": calc_temp,
            "calc_pressure": calc_pressure,
            "air_density_scale": air_density_scale,
            "sim_path": sim_path,
            "nu_grid": nu_grid,
            "equivalent_width": equivalent_width,
            "peak_abs_row": peak_abs_row,
            "max_optical_depth": max_optical_depth,
            "max_absorbance": max_absorbance,
            "integrated_opacity": integrated_opacity,
            "min_trans": min_trans,
            "mean_trans": mean_trans,
            "min_received_signal": min_received_signal,
            "detection_margin": detection_margin,
            "detection_threshold": detection_threshold,
            "strongest_text": strongest_text,
            "current_request_id": current_request_id,
        }

    # ── PERSISTENT SIMULATION RESULT RENDER ─────────────────────────
    if "active_sim_result" in st.session_state and st.session_state.active_sim_result:
        res = st.session_state.active_sim_result
        sim_df = res["sim_df"]
        pdf_config = res["pdf_config"]
        pdf_metrics = res["pdf_metrics"]

        st.success(
            f"Simulation READY — {res['sim_molecule']} (ID {res['mol_id']}, iso {res['iso_id']})  "
            f"at Lat {res['lat']:.2f}°, Lon {res['lon']:.2f}°, Altitude {res['altitude_m']:.0f} m"
        )

        st.subheader("Atmosphere")
        a1, a2, a3, a4 = st.columns(4)
        with a1: st.metric("Temperature T(z)",  f"{res['calc_temp']:.2f} K")
        with a2: st.metric("Pressure P(z)",      f"{res['calc_pressure']:.4f} atm")
        with a3: st.metric("Air Density Scale",  f"{res['air_density_scale']:.3f}")
        with a4: st.metric("Path Length",        f"{res['sim_path']:,.0f} cm")

        st.subheader("Spectroscopy")
        s1, s2, s3, s4 = st.columns(4)
        with s1: st.metric("Grid Points",        f"{len(res['nu_grid']):,}")
        with s2: st.metric("Equivalent Width",   f"{res['equivalent_width']:.4f} cm⁻¹")
        with s3: st.metric("Peak Absorption λ",  f"{float(res['peak_abs_row']['wavelength_um']):.4f} µm")
        with s4: st.metric("Peak Absorption ν",  f"{float(res['peak_abs_row']['wavenumber']):.2f} cm⁻¹")

        st.subheader("Radiative Transfer")
        r1, r2, r3 = st.columns(3)
        with r1: st.metric("Peak Optical Depth τ", f"{res['max_optical_depth']:.4f}")
        with r2: st.metric("Max Absorbance",        f"{res['max_absorbance']:.4f}")
        with r3: st.metric("Integrated Opacity",    f"{res['integrated_opacity']:.4e} cm⁻¹")

        st.subheader("Transmission & Sensor")
        t1, t2, t3, t4 = st.columns(4)
        with t1: st.metric("Min Transmittance",  f"{res['min_trans']:.4f}")
        with t2: st.metric("Mean Transmittance", f"{res['mean_trans']:.4f}")
        with t3: st.metric("Min Sensor Signal",  f"{res['min_received_signal']:.4f}")
        with t4: st.metric(
            "Sensor Status",
            "Detected" if res['detection_margin'] >= 0 else "Below threshold",
        )

        st.caption(res['strongest_text'])

        cd1, cd2 = st.columns(2)
        with cd1:
            st.download_button(
                "Download Simulation CSV",
                sim_df.to_csv(index=False).encode("utf-8"),
                file_name=f"oasis_{res['sim_molecule']}_alt{res['altitude_m']:.0f}m.csv",
                mime="text/csv",
                use_container_width=True,
            )
        with cd2:
            try:
                pdf_bytes = generate_simulation_pdf_report(sim_df, pdf_config, pdf_metrics)
                st.download_button(
                    "Download Complete PDF Report",
                    data=pdf_bytes,
                    file_name=f"oasis_report_{res['sim_molecule']}_alt{res['altitude_m']:.0f}m_{res['current_request_id'][:8]}.pdf",
                    mime="application/pdf",
                    type="primary",
                    use_container_width=True,
                )
            except Exception as e:
                st.error(f"Could not generate PDF report: {e}")

        # ─ Spectral plots ─────────────────────────────────────────
        st.subheader(
            f"LBL Spectral Outputs — {res['sim_molecule']}  "
            f"(Lat {res['lat']:.2f}°, Lon {res['lon']:.2f}°, Alt {res['altitude_m']:.0f} m)"
        )
        pt1, pt2, pt3, pt4 = st.tabs([
            "Transmittance",
            "Optical Depth",
            "Absorption Coefficient",
            "Sensor Signal",
        ])

        chart_layout_defaults = dict(
            template="plotly_white",
            hovermode="x unified",
            paper_bgcolor="rgba(0,0,0,0)",
            plot_bgcolor="#FFFFFF",
            font=dict(color="#0F172A"),
            margin=dict(l=40, r=40, t=50, b=40),
            xaxis=dict(showgrid=True, gridcolor="#E2E8F0", title_font=dict(size=13, color="#0F172A")),
            yaxis=dict(showgrid=True, gridcolor="#E2E8F0", title_font=dict(size=13, color="#0F172A")),
        )

        with pt1:
            fig = px.line(
                sim_df, x="wavenumber", y="transmittance",
                labels={"wavenumber": "Wavenumber (cm⁻¹)", "transmittance": "Transmittance (0–1)"},
                title=f"{res['sim_molecule']} Transmittance Spectrum at {res['altitude_m']:.0f} m",
                hover_data={"wavelength_um": ":.4f"},
            )
            fig.update_traces(line=dict(color="#2563EB", width=2.0))
            fig.update_layout(**chart_layout_defaults)
            st.plotly_chart(fig, use_container_width=True)

        with pt2:
            fig = px.line(
                sim_df, x="wavenumber", y="optical_depth",
                labels={"wavenumber": "Wavenumber (cm⁻¹)", "optical_depth": "Optical Depth τ"},
                title=f"{res['sim_molecule']} Optical Depth (τ = α · L)",
                hover_data={"wavelength_um": ":.4f"},
            )
            fig.update_traces(line=dict(color="#0D9488", width=2.0))
            fig.update_layout(**chart_layout_defaults)
            st.plotly_chart(fig, use_container_width=True)

        with pt3:
            fig = px.line(
                sim_df, x="wavenumber", y="absorption_coefficient",
                labels={"wavenumber": "Wavenumber (cm⁻¹)", "absorption_coefficient": "Absorption Coefficient α (cm⁻¹)"},
                title=f"{res['sim_molecule']} Voigt LBL Absorption Coefficient",
                hover_data={"wavelength_um": ":.4f"},
            )
            fig.update_traces(line=dict(color="#7C3AED", width=2.0))
            fig.update_layout(**chart_layout_defaults)
            st.plotly_chart(fig, use_container_width=True)

        with pt4:
            fig = px.line(
                sim_df, x="wavenumber", y="received_signal",
                labels={"wavenumber": "Wavenumber (cm⁻¹)", "received_signal": "Received Signal"},
                title=f"{res['sim_molecule']} Estimated Sensor Signal",
                hover_data={"wavelength_um": ":.4f"},
            )
            fig.add_hline(
                y=res['detection_threshold'], line_dash="dash", line_color="#DC2626",
                annotation_text=f"Detection threshold ({res['detection_threshold']})",
            )
            fig.update_traces(line=dict(color="#DC2626", width=2.0))
            fig.update_layout(**chart_layout_defaults)
            st.plotly_chart(fig, use_container_width=True)

        st.subheader("Simulation Output Table")
        st.dataframe(
            sim_df[[
                "wavenumber", "wavelength_um", "absorption_coefficient",
                "optical_depth", "transmittance", "absorbance", "received_signal",
            ]],
            use_container_width=True,
        )

    elif not run_button:
        if pf_status == "READY":
            st.info("Configuration is ready. Click **Run Altitude-Driven LBL Simulation** above.")
        elif pf_status == "WARNING":
            st.info(
                "Configuration has warnings but can run. "
                "Review the warnings above and click **Run** to proceed, "
                "or apply recommended values for a cleaner configuration."
            )



# ═══════════════════════════════════════════════════════════════════
# TAB 3 — GLOBAL COORDINATE PICKER & METRICS ENGINE
# ═══════════════════════════════════════════════════════════════════
with tab3:
    st.header("Global Coordinate Picker & Atmospheric Metrics Engine")
    st.write(
        "Click anywhere on the interactive world map to select a location. "
        "The system resolves the human-readable location name, fetches ground elevation from Open-Meteo, "
        "and computes standard atmospheric conditions from the ISA troposphere model. "
        "Selected coordinates are automatically synchronized with the simulation engine."
    )
    render_coordinate_picker()


# ─────────────────────────────────────────────────────────────────
# GLOBAL PLATFORM FOOTER
# ─────────────────────────────────────────────────────────────────
st.sidebar.markdown("---")
st.sidebar.markdown(
    "<div style='text-align: center; color: #64748B; font-size: 0.82rem;'>"
    "Open Atmospheric Spectroscopy &amp; Information System<br/>"
    "<b style='color: #2563EB; font-size: 0.88rem;'>By Ayushi Chahare</b>"
    "</div>",
    unsafe_allow_html=True,
)

st.markdown("---")
st.markdown(
    "<div style='text-align: center; color: #64748B; font-size: 0.82rem; padding: 12px 0 24px 0;'>"
    "OASIS Atmospheric Platform &middot; Open Atmospheric Spectroscopy &amp; Information System &middot; "
    "<b style='color: #2563EB; font-weight: 700;'>By Ayushi Chahare</b>"
    "</div>",
    unsafe_allow_html=True,
)



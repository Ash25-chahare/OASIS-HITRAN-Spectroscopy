# OASIS

OASIS is a Streamlit prototype for exploring local HITRAN/HAPI spectral line data, running altitude-driven line-by-line simulations, and inspecting station-level atmospheric metrics.

## What Is Included

- SQLite spectral database explorer
- HITRAN/HAPI line table cache in `data_cache`
- All 61 HITRAN molecule IDs selectable for LBL simulation
- On-demand HAPI fetch/cache for missing molecule and wavenumber ranges
- Molecule-specific LBL simulation with gas mole-fraction scaling
- Local validation summary before each simulation
- Wavenumber and wavelength input modes
- Molecule-specific suggested spectral windows for common atmospheric gases
- Atmosphere, spectroscopy, radiative-transfer, transmission, and sensor output sections
- Cursor-based global map picker with live latitude/longitude, click-to-select, elevation lookup, place lookup, favorites, and CSV export
- Interactive station reference map with hover labels and point selection
- Custom station/elevation metrics and map
- CSV exports for database and simulation results

## Run The App

From PowerShell:

```powershell
cd C:\Documents\OASIS_Project\backend
.\venv\Scripts\python.exe -m streamlit run app.py
```

This avoids PowerShell activation-policy issues.

## Rebuild The Local Database

```powershell
cd C:\Documents\OASIS_Project\backend
.\venv\Scripts\python.exe oasis.py
```

By default, the pipeline fetches isotopologue `1` for each configured molecule from `2200` to `2400 cm^-1`. The app can also fetch missing selected ranges on demand when you run a simulation.

You can change the fetch range without editing code:

```powershell
$env:OASIS_FETCH_START = "600"
$env:OASIS_FETCH_END = "800"
.\venv\Scripts\python.exe oasis.py
```

## Accuracy Notes

HITRAN line parameters come from the official HITRAN/HAPI data source. OASIS calculations are only as accurate as:

- the selected molecule and isotopologue
- the selected wavenumber range
- the local HAPI cache
- the temperature, pressure, and path length assumptions

Pressure and temperature are location/altitude values, so they can be identical across molecules. Absorption, transmittance, line count, line intensity, equivalent width, and optical depth should vary by molecule.

## LBL Model

The simulator uses HITRAN line parameters through HAPI. For the selected molecule, isotopologue, pressure, temperature, mole fraction, path length, and spectral window, HAPI calculates a Voigt line-shape absorption coefficient. OASIS then uses that coefficient to calculate transmittance over the path length.

Use narrow ranges for faster demos. Very wide ranges or many missing molecules may be slow because the app has to download and process new HAPI data.

Different molecules absorb strongly in different spectral windows. The simulator suggests common bands for major gases, then converts wavelength input back to wavenumber because HAPI performs calculations in `cm^-1`.

Dashboard flow:

```text
User spectral inputs + location inputs
-> altitude atmosphere model: T(z), P(z), air density
-> HITRAN/HAPI molecular line data
-> Voigt LBL spectroscopy model
-> radiative transfer: optical depth and absorbance
-> transmission spectrum
-> estimated sensor signal
```

## Global Station Map

The station engine has two map views.

The global cursor picker lets you move the cursor anywhere on the map to see latitude and longitude. Click a map point to lock that coordinate, fetch terrain elevation, estimate temperature and pressure, perform a place lookup, and save the point as a browser favorite.

The station reference map shows known/custom stations. Hover over a station to inspect latitude, longitude, altitude, estimated pressure, temperature, air-density scale, region, and tracked gases. Click a station point to select it in the inspector.

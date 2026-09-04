"""
OASIS Official Frontend Application Root
==========================================
Primary entrypoint for building, testing, and deploying the OASIS web user interface.
"""
import os
import sys
from pathlib import Path

# Add project root and backend directory to Python module search path
PROJECT_ROOT = Path(__file__).resolve().parent.parent
BACKEND_DIR = PROJECT_ROOT / "backend"

if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))
if str(BACKEND_DIR) not in sys.path:
    sys.path.insert(0, str(BACKEND_DIR))

# Execute backend application on every Streamlit script execution pass
import runpy
backend_app_path = BACKEND_DIR / "app.py"
runpy.run_path(str(backend_app_path), run_name="__main__")


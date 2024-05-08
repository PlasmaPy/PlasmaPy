"""Utility script for pulling ionization energy data from NIST and exporting it to a JSON file.

This script retrieves ionization energy data for elements from the NIST website,
formats it, and saves it as a JSON file for inclusion in the PlasmaPy package.
This script is provided for reference and is not intended to be part of the PlasmaPy package.
"""

import json
import logging
import time
from io import StringIO
from pathlib import Path

import pandas as pd
import requests

# Updated list of element symbols including Deuterium
elements = [
    "H",
    "D",
    "He",
    "Li",
    "Be",
    "B",
    "C",
    "N",
    "O",
    "F",
    "Ne",
    "Na",
    "Mg",
    "Al",
    "Si",
    "P",
    "S",
    "Cl",
    "Ar",
    "K",
    "Ca",
    "Sc",
    "Ti",
    "V",
    "Cr",
    "Mn",
    "Fe",
    "Co",
    "Ni",
    "Cu",
    "Zn",
    "Ga",
    "Ge",
    "As",
    "Se",
    "Br",
    "Kr",
    "Rb",
    "Sr",
    "Y",
    "Zr",
    "Nb",
    "Mo",
    "Tc",
    "Ru",
    "Rh",
    "Pd",
    "Ag",
    "Cd",
    "In",
    "Sn",
    "Sb",
    "Te",
    "I",
    "Xe",
    "Cs",
    "Ba",
    "La",
    "Ce",
    "Pr",
    "Nd",
    "Pm",
    "Sm",
    "Eu",
    "Gd",
    "Tb",
    "Dy",
    "Ho",
    "Er",
    "Tm",
    "Yb",
    "Lu",
    "Hf",
    "Ta",
    "W",
    "Re",
    "Os",
    "Ir",
    "Pt",
    "Au",
    "Hg",
    "Tl",
    "Pb",
    "Bi",
    "Po",
    "At",
    "Rn",
    "Fr",
    "Ra",
    "Ac",
    "Th",
    "Pa",
    "U",
    "Np",
    "Pu",
    "Am",
    "Cm",
    "Bk",
    "Cf",
    "Es",
    "Fm",
    "Md",
    "No",
    "Lr",
    "Rf",
    "Db",
    "Sg",
    "Bh",
    "Hs",
    "Mt",
    "Ds",
    "Rg",
    "Cn",
    "Nh",
    "Fl",
    "Mc",
    "Lv",
    "Ts",
    "Og",
]

# Dictionary to hold all ionization data
ionization_data = {}

# Base URL for the NIST ionization energy data
base_url = "https://physics.nist.gov/cgi-bin/ASD/ie.pl"


def add_to_dict(row) -> None:
    """Adds ionization data from a row of the DataFrame to the ionization_data dictionary."""
    ionization_data[row["ion"]] = {"ionization energy": row["ionization_energy"]}


# Iterate through each element
for element in elements:
    params = {
        "spectra": element,
        "submit": "Retrieve Data",
        "units": 1,
        "format": 2,
        "order": 0,
        "at_num_out": "on",
        "sp_name_out": "on",
        "ion_charge_out": "on",
        "el_name_out": "on",
        "seq_out": "on",
        "shells_out": "on",
        "level_out": "on",
        "ion_conf_out": "on",
        "e_out": 0,
        "unc_out": "on",
        "biblio": "on",
    }

    # Send GET request
    response = requests.get(base_url, params=params, timeout=10)

    try:
        text = response.text

        text = text.replace(",\n", "\n")

        # Use StringIO to handle CSV data as a string
        data = pd.read_csv(StringIO(response.text))

        # If there's a column named Ionization Energy (eV) (b), rename it to Ionization Energy (eV)
        if "Ionization Energy (eV) (b)" in data.columns:
            data = data.rename(
                columns={"Ionization Energy (eV) (b)": "Ionization Energy (eV)"}
            )

        # If there's a column named Ionization Energy (b) (eV), rename it to Ionization Energy (eV)
        if "Ionization Energy (b) (eV)" in data.columns:
            data = data.rename(
                columns={"Ionization Energy (b) (eV)": "Ionization Energy (eV)"}
            )

        # Drop all columns except for Sp. Name, Ion Charge, and Ionization Energy
        data = data[["Sp. Name", "Ion Charge", "Ionization Energy (eV)"]]

        # Remove quotes and = from all row values
        data["Sp. Name"] = data["Sp. Name"].str.replace('"', "")
        data["Sp. Name"] = data["Sp. Name"].str.replace("=", "")
        data["Ion Charge"] = data["Ion Charge"].str.replace('"', "")
        data["Ion Charge"] = data["Ion Charge"].str.replace("=", "")
        data["Ionization Energy (eV)"] = data["Ionization Energy (eV)"].str.replace(
            '"', ""
        )
        data["Ionization Energy (eV)"] = data["Ionization Energy (eV)"].str.replace(
            "=", ""
        )

        # In the name column, split the string by spaces and take the first element
        data["Sp. Name"] = data["Sp. Name"].str.split().str[0]

        # If the ion charge is not 0, append it to the name

        data["Sp. Name"] = (
            data["Sp. Name"] + " " + data["Ion Charge"].str.replace("+", "") + "+"
        )

        # If the ion charge is 0, remove the space and 0
        data["Sp. Name"] = data["Sp. Name"].str.replace(" 0+", "")

        # Rename Sp. Name to Ion
        data = data.rename(columns={"Sp. Name": "ion"})

        # Drop the ion charge column
        data = data.drop(columns=["Ion Charge"])

        # Try to convert Ionization Energy (eV) to float, and remove rows where it fails
        data["Ionization Energy (eV)"] = pd.to_numeric(
            data["Ionization Energy (eV)"], errors="coerce"
        )

        # Remove rows where Ionization Energy (eV) is NaN
        data = data.dropna(subset=["Ionization Energy (eV)"])

        # Convert Ionization Energy (eV) to float
        data["Ionization Energy (eV)"] = data["Ionization Energy (eV)"].astype(float)

        # Rename Ionization Energy (eV) to ionization energy
        data = data.rename(columns={"Ionization Energy (eV)": "ionization_energy"})

        # Add the data if ionization energy data is available; each ion is a separate record
        data.apply(add_to_dict, axis=1)

    except KeyError:
        logging.exception("Failed to parse data for %s", element)
    except requests.exceptions.RequestException:
        logging.exception("Failed to retrieve data for %s", element)

    # Delay of .5 seconds to avoid hitting rate limits or putting too much load on NIST
    time.sleep(0.5)

# Export the data to a JSON file
with Path.open(Path(__file__).parent / "ionization_energy.json", "w") as f:
    json.dump(ionization_data, f, indent=2)

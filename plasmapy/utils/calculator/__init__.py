""" 
Script and utilities to launch the plasma calculator
"""

import argparse
import os
import shlex
import subprocess

import plasmapy

def main():
    """
    Stub function for commandline tool that launches the plasma calculator notebook.
    """
    parser = argparse.ArgumentParser(description="Plasma calculator")
    parser.add_argument(
        "--port", type=int, default=8866, help="Port to run the notebook"
    )
    parser.add_argument(
        "--dark", action="store_true", help="Turn on dark mode, reduces eye strain"
    )
    parser.add_argument(
        "--no-browser", action="store_true", help="Do not open the browser"
    )

    module_path = plasmapy.__path__[0]
    computed_calculator_path = os.path.join(
        module_path, "utils", "calculator", "plasma_calculator.ipynb"
    )

    args = parser.parse_args()
    theme = "dark" if args.dark else "light"
    no_browser = "--no-browser" if args.no_browser else ""

    command = f"voila {no_browser} --port={args.port} --theme={theme} {computed_calculator_path} \
        --VoilaConfiguration.file_whitelist favicon.ico"
    try:
        subprocess.call(shlex.split(command))
    except KeyboardInterrupt:
        print("Stopping calculator! Bye")

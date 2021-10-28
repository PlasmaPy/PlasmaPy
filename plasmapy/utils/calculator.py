""" Script file to launch the notebook as a standalone clean app """

import argparse
import os
import shlex
import subprocess

import plasmapy


def main():
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
    calculator_path = os.path.dirname(module_path)
    computed_calculator_path = os.path.join(
        calculator_path, "docs", "notebooks", "plasma_calculator.ipynb"
    )

    args = parser.parse_args()
    theme = "dark" if args.dark else "light"
    no_browser = "--no-browser" if args.no_browser else ""

    command = f"voila {no_browser} --port={args.port} --theme={theme} {computed_calculator_path} \
        --VoilaConfiguration.file_whitelist=\"['favicon.ico']\""
    subprocess.Popen(shlex.split(command))

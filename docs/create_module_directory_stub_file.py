import glob
import os

from pathlib import Path


def main():
    directories = [str(d[0]).removeprefix("../") for d in os.walk("../plasmapy")]

    directories = [
        d.replace("/", ".")
        for d in directories
        if not d.endswith(("tests")) and "." not in d and "/_" not in d
    ]

    files = glob.glob("plasmapy/**/[a-zA-Z]*.py", root_dir="..")

    modules = sorted([*files, *directories])

    filepath = Path("./modules.rst")

    with filepath.open(mode="w") as f:
        preamble = ("Modules\n", "=======\n", "\n", ":orphan:\n", "\n")

        f.writelines(preamble)

        for file in modules:
            if "/_" in file:
                continue
            module = file.strip().removesuffix(".py").replace("/", ".")
            f.write(f"* `{module}`" + "\n")


if __name__ == "__main__":
    main()

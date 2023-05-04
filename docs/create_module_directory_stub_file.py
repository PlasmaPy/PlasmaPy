import glob
import os

from pathlib import Path


def get_subpackages(path: str) -> list[str]:
    directories = [str(d[0]).removeprefix("../") for d in os.walk(path)]

    return [
        d.replace("/", ".").strip()
        for d in directories
        if not d.endswith(("tests")) and "." not in d and "/_" not in d
    ]


def get_modules(path: str) -> list[str]:
    modules = list(glob.glob(f"{path}/**/[a-zA-Z]*.py"))

    return [
        m.removeprefix("../").removesuffix(".py").replace("/", ".")
        for m in modules
        if "tests/test_" not in m and "/_" not in m
    ]


def main():
    directories = [*get_subpackages("../plasmapy"), *get_subpackages("plasmapy_sphinx")]
    files = [*get_modules("../plasmapy"), *get_modules("plasmapy_sphinx")]
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

        f.writelines("\n")


if __name__ == "__main__":
    main()

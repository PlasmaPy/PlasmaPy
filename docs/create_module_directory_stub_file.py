import glob

from pathlib import Path


def main():
    files = glob.glob("plasmapy/**/[a-zA-Z]*.py", root_dir="..")
    filepath = Path("./modules.rst")

    with filepath.open(mode="w") as f:
        preamble = ("Modules\n", "=======\n", "\n", ":orphan:\n", "\n")

        f.writelines(preamble)

        for file in files:
            if "/_" in file:
                continue
            module = file.strip().removesuffix(".py").replace("/", ".")
            f.write(f"* `{module}`" + "\n")


if __name__ == "__main__":
    main()

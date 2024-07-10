"""Regenerate :file:`docs/changelog/index.rst`."""

import pathlib

repo = pathlib.Path(__file__).parent.parent
changelog_directory = repo / "changelog"
changelog_docs_directory = repo / "docs" / "changelog"
changelog_index_file = changelog_docs_directory / "index2.rst"

preamble = r""".. _changelog:

#########
Changelog
#########

The following pages describe the changes made during each release of
PlasmaPy.

.. toctree::
   :maxdepth: 1

"""


def there_are_unreleased_changes() -> bool:
    path = pathlib.Path(changelog_directory)
    return bool(path.glob("[1-9]*[0-9].*.rst"))


def get_release_notes_stems():
    files = changelog_docs_directory.glob("[0-9]*.[1-9]*.[0-9]*.rst")

    stems = [file.stem for file in files]

    def sort_key(stem: str) -> tuple[int, ...]:
        parts = stem.split(".")
        return tuple(int(part) for part in parts)

    return sorted(stems, key=sort_key, reverse=True)


def main():
    changelog_index_file.open(mode="w", encoding="utf-8")

    with changelog_index_file.open("w") as f:
        f.write(preamble)

        if there_are_unreleased_changes():
            f.write("   dev\n")

        stems = get_release_notes_stems()
        for stem in stems:
            f.write(f"   {stem}\n")

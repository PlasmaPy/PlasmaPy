"""
Regenerate the index of changelogs for releases.

Only include unreleased changes via `dev.rst` if there are changelog
entries corresponding to changes that have not been released.
"""

import pathlib

repo = pathlib.Path(__file__).parent.parent
changelog_directory = repo / "changelog"
changelog_docs_directory = repo / "docs" / "changelog"
changelog_index_file = changelog_docs_directory / "index.rst"

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
    """
    Return `True` if there are any news fragments in the top-level
    `changelog` directory corresponding to unreleased changes, and
    `False` otherwise.
    """
    path = pathlib.Path(changelog_directory)
    # path.glob() returns a generator, so make it into a list and then a bool
    return bool(list(path.glob("[1-9]*[0-9].*.rst")))


def get_release_notes_stems() -> list[str]:
    """
    Find the stems of filenames for the changelog files corresponding to
    releases.
    """
    files = changelog_docs_directory.glob("[0-9]*.[1-9]*.[0-9]*.rst")

    stems = [file.stem for file in files]

    def sort_key(stem: str) -> tuple[int, ...]:
        parts = stem.split(".")
        return tuple(int(part) for part in parts)

    return sorted(stems, key=sort_key, reverse=True)


def main() -> None:
    """
    Regenerate the index file for release changelog files. Include
    unreleased changes only when there are unreleased changes.
    """
    with changelog_index_file.open(mode="w", encoding="utf-8") as f:
        f.write(preamble)

        if there_are_unreleased_changes():
            f.write("   dev\n")

        stems = get_release_notes_stems()
        for stem in stems:
            f.write(f"   {stem}\n")

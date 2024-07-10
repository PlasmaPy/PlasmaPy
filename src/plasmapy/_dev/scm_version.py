# Try to use setuptools_scm to get the current version; this is only used
# in development installations from the git repository.
import os.path as pth
from datetime import date

try:
    from setuptools_scm import get_version

    version = get_version(
        root=pth.join("..", "..", ".."),
        relative_to=__file__,
    )

    # get_version infers the next version from the most recent tag on
    # `main`. However, we have been tagging versions on the release
    # branches rather than on `main`. We extract the git hash from
    # get_version, and reconstruct a date-based dev version. There's
    # undoubtedly a more elegant way to do this, but the docstring for
    # get_version doesn't describe how to use it.

    git_hash_and_date = version.split("+")[-1]
    git_hash = git_hash_and_date.split(".")[0]
    today = date.today()  # noqa: DTZ011
    version = f"{today.year}.{today.month}.{today.day}.dev+{git_hash}"

except Exception as e:
    raise ImportError("Unable to get version using setuptools_scm.") from e

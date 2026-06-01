__all__ = ["decorator"]

from collections.abc import Callable
from typing import Any, TypeVar

_F = TypeVar("_F", bound=Callable[..., Any])

def decorator(  # noqa: UP047
    wrapper: _F,
    enabled: Any = None,  # noqa: ANN401
    adapter: Any = None,  # noqa: ANN401
    proxy: Any = None,  # noqa: ANN401
) -> _F: ...

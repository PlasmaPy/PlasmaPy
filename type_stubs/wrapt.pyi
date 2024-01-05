__all__ = ["decorator"]

from typing import Any, Callable, TypeVar

_F = TypeVar("_F", bound=Callable[..., Any])

def decorator(
    wrapper: _F,
    enabled: Any = None,
    adapter: Any = None,
    proxy: Any = None,
) -> _F: ...

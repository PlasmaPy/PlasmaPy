"""
Generic classes for decorators that process arguments based on
annotations.
"""

__all__ = ["GenericDecorator"]

from abc import ABC, abstractmethod


class GenericDecorator(ABC):
    def __init__(self):
        pass

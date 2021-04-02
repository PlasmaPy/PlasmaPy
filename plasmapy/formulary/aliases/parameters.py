import importlib
import inspect
import sys

module = importlib.import_module(__name__.replace(".aliases", ""))
__all__ = []
for name, function in inspect.getmembers(module):
    if not name.startswith("_") and name.endswith("_"):
        # print(f"{name=} fits!")
        __all__.append(name)
        setattr(sys.modules[__name__], name, function)
    # else:
    #     print(f"\t{name=} SKIPPED!")
__all__.sort()
# print(__all__)

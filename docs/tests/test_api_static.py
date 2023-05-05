import os


def _is_valid_dir(path: str) -> bool:
    """
    Check if a directory is valid based on if it contains .py files and
    if it does not start with '_'.
    """
    if any(file.startswith("test_") for file in os.listdir(path)):
        return False
    py_files = [
        file
        for file in os.listdir(path)
        if file.endswith(".py")
        and not file.startswith("_")
        and file not in ("conftest.py", "test_conftest.py")
    ]
    return bool(py_files)


def _get_py_files(root):
    """Return a list of .py files in a directory"""
    return [
        file
        for file in os.listdir(root)
        if file.endswith(".py")
        and not file.startswith("_")
        and file not in ("conftest.py", "test_conftest.py")
    ]


def _get_py_namespaces(path):
    namespaces = []
    for root, dirs, files in os.walk(path):
        # Exclude directories that begin with an underscore
        dirs[:] = [d for d in dirs if not d.startswith("_")]
        if _is_valid_dir(root):
            rel_path = os.path.relpath(root, path)
            namespace = (
                f"{os.path.basename(path)}.{rel_path.replace('/', '.')}"
                if rel_path != "."
                else os.path.basename(path)
            )
            namespaces.append(namespace)
            for file in _get_py_files(root):
                rel_path = os.path.relpath(os.path.join(root, file), path)
                namespace = f"{os.path.basename(path)}.{os.path.splitext(rel_path)[0].replace('/', '.')}"
                namespaces.append(namespace)
    return sorted(namespaces)

from setuptools import setup
from setuptools_rust import Binding, RustExtension

setup(
    rust_extensions=[
        RustExtension(
            "plasmapy.rust",
            path="Cargo.toml",
            binding=Binding.PyO3
        )
    ],
)

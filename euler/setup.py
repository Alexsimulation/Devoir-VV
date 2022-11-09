from glob import glob
from setuptools import setup, Extension
from pybind11.setup_helpers import Pybind11Extension


# Run with setup.py build --force

__version__ = "0.0.1"


ext_modules = [
    Pybind11Extension(
        "euler",
        ["./src/python.cpp"],  # Sort source files for reproducibility
        include_dirs=["./"],
        extra_compile_args = ["-std=c++17", "-fopenmp", "-fopenmp-simd", "-Ofast", "-w", "-DEULER_PYTHON_MODULE"],
        extra_link_args=['-lgomp']
    )
]

setup(
    name="euler",
    version=__version__,
    author="Alexis Angers",
    author_email="alexis.angers@polymtl.ca",
    url="https://github.com/AER8875-2022/Euler-Dev/tree/master/Alexis/euler",
    description="Python API for Euler equations solver",
    long_description="",
    ext_modules=ext_modules
)


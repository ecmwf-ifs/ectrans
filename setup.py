import os
import ast
from skbuild import setup

_version_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), "VERSION")
with open(_version_file, "r") as f:
    __version__ = f.read().strip()

setup(
    name="ectrans4py",
    version=__version__,
    packages=['ectrans4py'],
    cmake_minimum_required_version="3.13",
    cmake_args=[
        '-DENABLE_ETRANS=ON',
        '-DENABLE_ECTRANS4PY=ON',
        '-DENABLE_SINGLE_PRECISION=OFF',
        '-DENABLE_OMP=OFF',
    ],
    package_dir={"": "src"},
    cmake_install_dir="src/ectrans4py",
    setup_requires=["scikit-build", "setuptools"],
    install_requires=["numpy", "ctypesforfortran==1.1.3"],
)

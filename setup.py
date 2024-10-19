import os
import re
import ast
from skbuild import setup

def get_version():   # remove this part 
    version_file = os.path.join("src", "ectrans4py", "__init__.py")
    with open(version_file, "r", encoding="utf-8") as f:
        content = f.read()
        version_match = re.search(r"^__version__\s*=\s*['\"]([^'\"]*)['\"]", content, re.M)
        if version_match:
            return version_match.group(1)
        raise RuntimeError("Unable to find version string.")

version=get_version()
# ectrans4py package : 
setup(
    name="ectrans4py",
    version=version, 
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

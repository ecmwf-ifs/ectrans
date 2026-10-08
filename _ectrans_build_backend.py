# (C) Copyright 2026- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

"""In-tree PEP 517 build backend for ectrans4py.

This is a thin wrapper around ``scikit_build_core.build``. Before a wheel is
built it runs ``ectrans-bundle create`` so the ecbuild/fiat dependencies are
cloned into a throwaway directory (reusing the existing ecbundle machinery
rather than duplicating the git URLs/versions). scikit-build-core is then
pointed at the generated bundle superbuild by overriding ``cmake.source-dir``.

The bundle is created in a temporary directory rather than
``package/bundle/source`` so that a bundle the user may already be working on in
the checkout is never clobbered, and the clone is discarded once the wheel has
been built.
"""

import os
import shutil
import subprocess
import tempfile
from contextlib import contextmanager
from pathlib import Path

from scikit_build_core import build as _skbuild

_ROOT = Path(__file__).parent.resolve()
_BUNDLE_DIR = _ROOT / "package" / "bundle"


@contextmanager
def _bundle_sources():
    """Create the ecbuild/fiat/ectrans superbuild in a temporary directory.

    Yields the source directory to point ``cmake.source-dir`` at, and removes
    the whole checkout on exit.
    """
    tmp = Path(tempfile.mkdtemp(prefix="ectrans4py-bundle-"))
    src = tmp / "source"
    try:
        # cwd=_BUNDLE_DIR so ecbundle discovers bundle.yml there.
        subprocess.run(
            ["./ectrans-bundle", "create", "--src-dir", str(src)],
            cwd=_BUNDLE_DIR,
            check=True,
        )

        # bundle.yml uses `dir: $PWD` for ectrans, but $PWD resolves to the
        # wrapper script's directory when create runs here, not this repo. Point
        # the symlink at the actual checkout so add_subdirectory(ectrans) works.
        ectrans_link = src / "ectrans"
        if ectrans_link.is_symlink() or ectrans_link.exists():
            ectrans_link.unlink()
        os.symlink(_ROOT, ectrans_link)

        yield src
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def _with_source_dir(config_settings, src):
    settings = dict(config_settings or {})
    settings["cmake.source-dir"] = str(src)
    return settings


def build_wheel(wheel_directory, config_settings=None, metadata_directory=None):
    with _bundle_sources() as src:
        return _skbuild.build_wheel(
            wheel_directory, _with_source_dir(config_settings, src), metadata_directory
        )


def build_editable(wheel_directory, config_settings=None, metadata_directory=None):
    with _bundle_sources() as src:
        return _skbuild.build_editable(
            wheel_directory, _with_source_dir(config_settings, src), metadata_directory
        )


def build_sdist(sdist_directory, config_settings=None):
    return _skbuild.build_sdist(sdist_directory, config_settings)


# Remaining PEP 517 hooks are implemented by scikit-build-core and delegated as-is.
get_requires_for_build_wheel = _skbuild.get_requires_for_build_wheel
get_requires_for_build_sdist = _skbuild.get_requires_for_build_sdist
get_requires_for_build_editable = _skbuild.get_requires_for_build_editable
prepare_metadata_for_build_wheel = _skbuild.prepare_metadata_for_build_wheel
prepare_metadata_for_build_editable = _skbuild.prepare_metadata_for_build_editable

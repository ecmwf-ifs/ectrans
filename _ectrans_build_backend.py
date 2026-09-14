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
cloned into ``package/bundle/source`` (reusing the existing ecbundle machinery
rather than duplicating the git URLs/versions). scikit-build-core is then
pointed at the generated bundle superbuild via ``cmake.source-dir``.
"""

import subprocess
from pathlib import Path

from scikit_build_core import build as _skbuild

_ROOT = Path(__file__).parent.resolve()
_BUNDLE_DIR = _ROOT / "package" / "bundle"
_BUNDLE_SOURCE = _BUNDLE_DIR / "source"


def _ensure_bundle_sources():
    """Clone ecbuild/fiat via ectrans-bundle if they are not already present."""
    if (_BUNDLE_SOURCE / "ecbuild").is_dir() and (_BUNDLE_SOURCE / "fiat").is_dir():
        return
    subprocess.run(
        ["./ectrans-bundle", "create"],
        cwd=_BUNDLE_DIR,
        check=True,
    )


def build_wheel(wheel_directory, config_settings=None, metadata_directory=None):
    _ensure_bundle_sources()
    return _skbuild.build_wheel(wheel_directory, config_settings, metadata_directory)


def build_editable(wheel_directory, config_settings=None, metadata_directory=None):
    _ensure_bundle_sources()
    return _skbuild.build_editable(wheel_directory, config_settings, metadata_directory)


def build_sdist(sdist_directory, config_settings=None):
    return _skbuild.build_sdist(sdist_directory, config_settings)


# Remaining PEP 517 hooks are implemented by scikit-build-core and delegated as-is.
get_requires_for_build_wheel = _skbuild.get_requires_for_build_wheel
get_requires_for_build_sdist = _skbuild.get_requires_for_build_sdist
get_requires_for_build_editable = _skbuild.get_requires_for_build_editable
prepare_metadata_for_build_wheel = _skbuild.prepare_metadata_for_build_wheel
prepare_metadata_for_build_editable = _skbuild.prepare_metadata_for_build_editable

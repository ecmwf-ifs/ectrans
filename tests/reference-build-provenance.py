#!/usr/bin/env python3

# (C) Copyright 2026- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

"""Collect build provenance for generated checksum references.

Intended usage
--------------
This script is designed to run as part of the reference-update flow, not as a
general-purpose system inventory tool. The expected caller is
`tests/ectrans-update-references.sh` (generated from
`tests/ectrans-update-references.sh.in`), which is itself invoked by
`.github/tools/reference-artifact-package.sh` when assembling the reference
artifact consumed by `.github/workflows/update-references.yml`.

Inputs
------
- `--build-dir`: exact CMake build directory for the ecTrans build. The
    directory must contain `CMakeCache.txt`.
- `--output`: JSON file path to write.
- `--anonymous` (optional): replace local ecbuild, fiat, and field_api roots
    in dependency metadata with their corresponding CMake variable names.
- `--dependency-ref NAME=REF` (optional, repeatable): requested refs for
    dependency version annotation.

Output
------
Writes one JSON document containing system/compiler/cmake/dependency metadata.
Before writing, string data is sanitized to reduce leakage of secret-like env
values.
"""

import argparse
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile


# These environment variables are considered sensitive for reproducibility and will be included in
# the provenance output if set.

ENVIRONMENT_VARIABLES = (
    "MKL_CBWR", # Controls MKL behavior for reproducibility
    "LIBSCI_ARCH_OVERRIDE", # Cray LibSci architecture override
)

def read_cmake_cache(build_dir):
    cache = {}
    cache_file = build_dir / "CMakeCache.txt"
    if not cache_file.is_file():
        return cache

    for line in cache_file.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line or line.startswith("//") or line.startswith("#") or "=" not in line:
            continue
        key_type, value = line.split("=", 1)
        if ":" not in key_type:
            continue
        key, _ = key_type.split(":", 1)
        cache[key] = value
    return cache


def read_manifest(build_dir):
    values = {}
    manifest = build_dir / "ectrans-reference-provenance.txt"
    if not manifest.is_file():
        return values

    for line in manifest.read_text(encoding="utf-8", errors="replace").splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            values[key] = value
    return values


def cmake_list(value):
    return [item for item in value.split(";") if item and not item.endswith("-NOTFOUND")]


def absolute_path(path):
    if not path:
        return None
    candidate = Path(path).expanduser()
    if not candidate.is_absolute():
        candidate = Path(shutil.which(str(candidate)) or candidate)
    return str(candidate.resolve())


def absolute_paths(paths):
    return sorted(set(absolute_path(path) for path in paths if path))


def cache_library_paths(cache, prefixes):
    paths = []
    for key, value in cache.items():
        if not key.startswith(prefixes) or not key.endswith(("_LIB", "_LIBRARY")):
            continue
        if value and not value.endswith("-NOTFOUND"):
            paths.extend(cmake_list(value))
    return absolute_paths(paths)


def installed_library_paths(root):
    if not root:
        return []
    root = Path(absolute_path(root))
    paths = []
    for directory in (root / "lib", root / "lib64"):
        if directory.is_dir():
            paths.extend(str(path) for path in directory.glob("lib*.*") if path.is_file())
    return absolute_paths(paths)


def library_root(library_paths):
    if not library_paths:
        return None
    library = Path(library_paths[0])
    for directory in (library.parent, *library.parents):
        if directory.name in ("lib", "lib64"):
            return absolute_path(directory.parent)
    return absolute_path(library.parent)


def first_value(*values):
    return next((value for value in values if value), None)


# Python 3.6 does not support subprocess text=True, so decode bytes explicitly.
def decode_output(output):
    return output.decode("utf-8", errors="replace").strip()


def command_output(command, cwd=None):
    try:
        return decode_output(subprocess.check_output(command, cwd=cwd, stderr=subprocess.DEVNULL))
    except (OSError, subprocess.CalledProcessError):
        return None


def sensitive_env_values():
    values = []
    for key, value in os.environ.items():
        if not value or len(value) < 4:
            continue
        if re.search(r"TOKEN|SECRET|PASSWORD|PASS|KEY", key, re.IGNORECASE):
            values.append(value)
    return sorted(set(values), key=len, reverse=True)


def redact_string(value, secrets):
    result = value

    for secret in secrets:
        result = result.replace(secret, "***")

    return result


def redact(value, secrets):
    if isinstance(value, dict):
        return {key: redact(item, secrets) for key, item in value.items()}
    if isinstance(value, list):
        return [redact(item, secrets) for item in value]
    if isinstance(value, str):
        return redact_string(value, secrets)
    return value


def anonymous_path(path, root, placeholder):
    if not path or not root:
        return path
    normalized_path = os.path.normpath(path)
    normalized_root = os.path.normpath(absolute_path(root))
    if normalized_path == normalized_root:
        return placeholder
    prefix = normalized_root + os.sep
    if normalized_path.startswith(prefix):
        return placeholder + os.sep + normalized_path[len(prefix):]
    return path


def anonymize_dependency(dependency_info, root, placeholder):
    if not root:
        return
    dependency_info["root"] = placeholder
    dependency_info["library_paths"] = [
        anonymous_path(path, root, placeholder) for path in dependency_info["library_paths"]
    ]


def git_release_tag(path):
    if not path:
        return None

    path = Path(path)
    for candidate in (path, *path.parents):
        if (candidate / ".git").exists():
            tags = command_output(["git", "tag", "--points-at", "HEAD"], candidate)
            if tags:
                release_tags = sorted(tag for tag in tags.splitlines() if re.fullmatch(r"v?\d+(?:\.\d+)+(?:[-+][0-9A-Za-z.-]+)?", tag))
                if release_tags:
                    return release_tags[-1]
            return None
    return None


def cache_dependency_version(cache, name):
    name_upper = name.upper()
    name_candidates = (name, name_upper)

    scalar_keys = []
    for candidate in name_candidates:
        scalar_keys.extend(
            (
                f"{candidate}_VERSION_STR",
                f"{candidate}_VERSION_STRING",
                f"{candidate}_VERSION",
            )
        )

    for key in scalar_keys:
        value = cache.get(key)
        if value and not value.endswith("-NOTFOUND"):
            return value

    for candidate in name_candidates:
        major = cache.get(f"{candidate}_VERSION_MAJOR")
        minor = cache.get(f"{candidate}_VERSION_MINOR")
        patch = cache.get(f"{candidate}_VERSION_PATCH")
        if all((major, minor, patch)):
            return f"{major}.{minor}.{patch}"

    return None


def manifest_dependency_version(manifest, name):
    name_upper = name.upper()
    keys = (
        f"{name_upper}_VERSION_STR",
        f"{name_upper}_VERSION_STRING",
        f"{name_upper}_VERSION",
    )
    for key in keys:
        value = manifest.get(key)
        if value and not value.endswith("-NOTFOUND"):
            return value
    return None


def dependency_version(name, dependency_refs, cache, manifest):
    repository = cache.get(f"{name.upper()}_SOURCE_DIR") or manifest.get(f"{name.upper()}_SOURCE_DIR") or cache.get(f"{name.upper()}_ROOT") or manifest.get(f"{name.upper()}_ROOT")
    return first_value(
        dependency_ref(name, dependency_refs),
        git_release_tag(repository),
        manifest_dependency_version(manifest, name),
        cache_dependency_version(cache, name),
    )


def dependency_ref(name, cli_refs):
    cli_value = cli_refs.get(name)
    if cli_value:
        return cli_value

    env_key = f"ECTRANS_DEP_REF_{re.sub(r'[^A-Za-z0-9]', '_', name).upper()}"
    env_value = os.environ.get(env_key)
    return env_value if env_value else None


def version_from_paths(paths, package_name):
    pattern = re.compile(rf"(?:^|[/_-]){re.escape(package_name)}[/_-](\d+(?:\.\d+)+(?:[-+][0-9A-Za-z.-]+)?)(?:/|$)", re.IGNORECASE)
    for path in paths:
        match = pattern.search(path)
        if match:
            return match.group(1)
    return None


def environment_info():
    return {name: os.environ[name] for name in ENVIRONMENT_VARIABLES if os.environ.get(name)}


def intel_mkl_version(root, library_paths):
    candidates = [root, os.environ.get("MKLROOT")]
    candidates.extend(library_paths)
    candidates.extend(module for module in os.environ.get("LOADEDMODULES", "").split(":") if "mkl" in module.lower())

    patterns = (
        r"compilers_and_libraries_(\d+(?:\.\d+)+)",
        r"parallel_studio_xe_(\d+_update\d+)",
        r"(?:^|[/_-])(?:intel[-_/])?mkl[/_-](\d+(?:\.\d+)+(?:[-+][0-9A-Za-z.-]+)?)(?:/|$)",
        r"(?:^|[/_-])oneapi[/_-]mkl[/_-](\d+(?:\.\d+)+(?:[-+][0-9A-Za-z.-]+)?)(?:/|$)",
    )

    for candidate in candidates:
        if not candidate:
            continue
        for pattern in patterns:
            match = re.search(pattern, candidate, re.IGNORECASE)
            if match:
                return match.group(1)
    return None


def compiler_info(cache, manifest, language):
    cmake_language = {"c": "C", "cxx": "CXX", "fortran": "Fortran", "cuda": "CUDA", "hip": "HIP"}[language]
    compiler = cache.get(f"CMAKE_{cmake_language}_COMPILER")
    if not compiler or compiler.endswith("-NOTFOUND"):
        return None
    compiler = absolute_path(compiler)

    build_type = cache.get("CMAKE_BUILD_TYPE", "").upper()
    version_output = command_output([compiler, "--version"])
    version_line = version_output.splitlines()[0] if version_output else ""
    compiler_id = first_value(cache.get(f"CMAKE_{cmake_language}_COMPILER_ID"), manifest.get(f"CMAKE_{cmake_language}_COMPILER_ID"))
    compiler_version = first_value(cache.get(f"CMAKE_{cmake_language}_COMPILER_VERSION"), manifest.get(f"CMAKE_{cmake_language}_COMPILER_VERSION"))
    if not compiler_id:
        if "NVHPC" in version_line or "nvfortran" in version_line:
            compiler_id = "NVHPC"
        elif "Intel" in version_line:
            compiler_id = "Intel"
        elif "clang" in version_line.lower():
            compiler_id = "Clang"
        elif "GNU" in version_line or "gcc" in version_line.lower() or "gfortran" in version_line.lower():
            compiler_id = "GNU"
    if not compiler_version:
        match = re.search(r"\b\d+(?:\.\d+)+(?:\.\d+)?\b", version_line)
        compiler_version = match.group(0) if match else None

    return {
        "path": compiler,
        "id": compiler_id,
        "version": compiler_version,
        "flags": cache.get(f"CMAKE_{cmake_language}_FLAGS", ""),
        "build_type_flags": cache.get(f"CMAKE_{cmake_language}_FLAGS_{build_type}", ""),
    }


def cmake_info(cache, manifest):
    version = manifest.get("CMAKE_VERSION")
    if not version:
        version_parts = [cache.get(f"CMAKE_CACHE_{part}_VERSION") for part in ("MAJOR", "MINOR", "PATCH")]
        version = ".".join(version_parts) if all(version_parts) else None
    if not version:
        version_output = command_output(["cmake", "--version"])
        match = re.search(r"cmake version\s+(\S+)", version_output or "")
        version = match.group(1) if match else None
    return {"version": version}


def lapack_version(compiler, library_paths):
    mkl_version = intel_mkl_version(library_root(library_paths), library_paths)
    if mkl_version:
        return mkl_version

    if not compiler or not library_paths:
        return None

    source = """program lapack_version
  implicit none
  integer :: major, minor, patch
  call ilaver(major, minor, patch)
  write (*, '(I0,".",I0,".",I0)') major, minor, patch
end program lapack_version
"""
    with tempfile.TemporaryDirectory() as temp_dir:
        source_path = Path(temp_dir) / "lapack_version.F90"
        executable = Path(temp_dir) / "lapack_version"
        source_path.write_text(source, encoding="utf-8")
        try:
            subprocess.run(
                [compiler, str(source_path), "-o", str(executable), *library_paths],
                check=True,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )
            return decode_output(subprocess.check_output([str(executable)], stderr=subprocess.DEVNULL))
        except (OSError, subprocess.CalledProcessError):
            return None


def system_info():
    cpu_info = None
    if sys.platform == "darwin":
        cpu_info = command_output(["system_profiler", "SPHardwareDataType", "-detailLevel", "mini"])
        cpu_brand = command_output(["sysctl", "-n", "machdep.cpu.brand_string"])
        if cpu_brand:
            cpu_info = "\n".join(filter(None, (cpu_info, f"machdep.cpu.brand_string: {cpu_brand}")))
    elif Path("/proc/cpuinfo").is_file():
        cpu_info = Path("/proc/cpuinfo").read_text(encoding="utf-8", errors="replace").split("\n\n", 1)[0]

    processor = None
    if cpu_info:
        match = re.search(r"^\s*(?:model name|Processor|Chip|machdep\.cpu\.brand_string)\s*[:=]\s*(.+)$", cpu_info, re.MULTILINE)
        if match:
            processor = match.group(1).strip()

    runner_os = os.environ.get("RUNNER_OS") or {"darwin": "macos"}.get(sys.platform, sys.platform)
    runner_arch = os.environ.get("RUNNER_ARCH") or os.uname().machine

    return {
        "os": runner_os,
        "arch": runner_arch,
        "processor": processor,
        "loaded_modules": [module for module in os.environ.get("LOADEDMODULES", "").split(":") if module],
    }


def dependency(name, vendor, version, root, library_paths, **extra):
    result = {
        "vendor": vendor,
        "version": version,
        "root": absolute_path(root),
        "library_paths": absolute_paths(library_paths),
    }
    result.update(extra)
    return result


def module_version(prefix):
    for module in os.environ.get("LOADEDMODULES", "").split(":"):
        if module == prefix or module.startswith(f"{prefix}/"):
            return module.partition("/")[2] or None
    return None


def lapack_vendor(library_paths, cache=None, manifest_libraries=None):
    if os.environ.get("ECTRANS_LAPACK_VENDOR"):
        return os.environ["ECTRANS_LAPACK_VENDOR"]
    if os.environ.get("MKLROOT"):
        return "Intel MKL"
    if os.environ.get("CRAY_LIBSCI_DIR") or os.environ.get("CRAY_LIBSCI_BASE_DIR"):
        return "Cray LibSci"

    framework_entries = [entry.lower() for entry in cmake_list(manifest_libraries or "")]
    if "-framework" in framework_entries and any("accelerate" == entry or "accelerate.framework" in entry for entry in framework_entries):
        return "Apple Accelerate"

    if os.environ.get("OPENBLAS_ROOT") or os.environ.get("OpenBLAS_DIR"):
        return "OpenBLAS"

    if cache and cache.get("BLA_VENDOR", "").lower() in ("apple", "accelerate"):
        return "Apple Accelerate"

    for path in library_paths:
        try:
            resolved = str(Path(path).expanduser().resolve()).lower()
        except OSError:
            resolved = str(Path(path).expanduser()).lower()
        name = Path(path).name.lower()
        if "openblas" in name or "openblas" in resolved:
            return "OpenBLAS"
        if "accelerate.framework" in name or "accelerate.framework" in resolved:
            return "Apple Accelerate"
        if "veclib.framework" in name or "veclib.framework" in resolved:
            return "Apple Accelerate"

    return "Netlib"


def fftw_vendor(library_paths, cache=None):
    if os.environ.get("ECTRANS_FFTW_VENDOR"):
        return os.environ["ECTRANS_FFTW_VENDOR"]
    if cache and cache.get("FFTW_ENABLE_MKL", "").upper() in ("ON", "TRUE", "YES", "1"):
        return "Intel MKL"
    if os.environ.get("MKLROOT"):
        return "Intel MKL"
    if os.environ.get("NVPL_ROOT") or os.environ.get("NVPL_INCLUDE_DIR") or os.environ.get("NVPL_INCLUDE_DIR_FFTW"):
        return "NVIDIA NVPL"
    if os.environ.get("CRAY_FFTW_PREFIX") or module_version("cray-fftw"):
        return "Cray FFTW"

    for path in library_paths:
        try:
            resolved = str(Path(path).expanduser().resolve()).lower()
        except OSError:
            resolved = str(Path(path).expanduser()).lower()
        name = Path(path).name.lower()
        text = f"{name} {resolved}"
        if "mkl" in text:
            return "Intel MKL"
        if "nvpl" in text:
            return "NVIDIA NVPL"
        if "cray" in text or "libsci" in text:
            return "Cray FFTW"
        if "fftw" in text:
            return "FFTW"

    return "FFTW-compatible"


def mpi_vendor(root, version_output):
    text = " ".join(value for value in (root, version_output) if value).lower()
    if "open" in text and "mpi" in text:
        return "Open MPI"
    if "mpich" in text:
        return "MPICH"
    if "intel" in text or "impi" in text:
        return "Intel MPI"
    if "cray" in text or "mpich" in os.environ.get("PE_ENV", "").lower():
        return "Cray MPICH"
    return "MPI"


def mpi_version(cache):
    executable = first_value(cache.get("MPIEXEC_EXECUTABLE"), cache.get("MPIEXEC"), cache.get("MPI_C_COMPILER"))
    output = command_output([executable, "--version"]) if executable else None
    if not output and cache.get("MPI_C_COMPILER"):
        output = command_output([cache["MPI_C_COMPILER"], "--version"])
    match = re.search(r"\b\d+(?:\.\d+)+(?:[-+][0-9A-Za-z.-]+)?\b", output or "")
    return (match.group(0) if match else None), output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build-dir", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--anonymous", action="store_true")
    parser.add_argument("--dependency-ref", action="append", default=[], metavar="NAME=REF")
    args = parser.parse_args()

    dependency_refs = {}
    for entry in args.dependency_ref:
        name, separator, ref = entry.partition("=")
        if not separator:
            parser.error(f"Invalid dependency reference: {entry}")
        dependency_refs[name] = ref

    try:
        build_dir = args.build_dir.expanduser().resolve()
    except OSError:
        build_dir = args.build_dir.expanduser().absolute()
    if not build_dir.is_dir():
        parser.error(f"Build directory does not exist: {build_dir}")
    if not (build_dir / "CMakeCache.txt").is_file():
        parser.error(f"Build directory is missing CMakeCache.txt: {build_dir}")

    cache = read_cmake_cache(build_dir)
    manifest = read_manifest(build_dir)
    compilers = {}
    for name in ("c", "cxx", "fortran", "cuda", "hip"):
        info = compiler_info(cache, manifest, name)
        if info:
            compilers[name] = info

    lapack_library_paths = cmake_list(manifest.get("LAPACK_LIBRARIES", ""))
    lapack_library_paths = [path for path in lapack_library_paths if Path(path).is_file()] or cache_library_paths(cache, ("LAPACK_", "BLAS_"))
    lapack_library_paths = absolute_paths(lapack_library_paths)
    fftw_library_paths = cmake_list(manifest.get("FFTW_LIBRARIES", ""))
    fftw_library_paths = [path for path in fftw_library_paths if Path(path).is_file()] or cache_library_paths(cache, ("FFTW_", "pkgcfg_lib_PKG_FFTW_"))
    fftw_library_paths = absolute_paths(fftw_library_paths)
    mpi_library_paths = cache_library_paths(cache, ("MPI_",))
    mpi_library_paths = absolute_paths(mpi_library_paths)

    fftw_root = first_value(library_root(fftw_library_paths), os.environ.get("FFTW_ROOT"), cache.get("FFTW_ROOT"))
    lapack_root = first_value(library_root(lapack_library_paths), cache.get("LAPACK_ROOT"))
    mpi_root = first_value(library_root(mpi_library_paths), os.environ.get("MPI_HOME"))
    fiat_root = first_value(cache.get("fiat_ROOT"), cache.get("fiat_DIR"))
    field_api_root = first_value(cache.get("field_api_ROOT"), cache.get("field_api_DIR"))
    ecbuild_root = cache.get("ecbuild_DIR")
    fortran_compiler = compilers.get("fortran", {}).get("path")
    lapack_dependency_vendor = lapack_vendor(lapack_library_paths, cache=cache, manifest_libraries=manifest.get("LAPACK_LIBRARIES", ""))
    fftw_dependency_vendor = fftw_vendor(fftw_library_paths, cache=cache)
    fftw_version = first_value(
        intel_mkl_version(fftw_root, fftw_library_paths) if fftw_dependency_vendor == "Intel MKL" else None,
        manifest.get("FFTW_VERSION"),
        cache.get("FFTW_VERSION"),
        manifest_dependency_version(manifest, "fftw"),
        version_from_paths(fftw_library_paths, "fftw"),
        module_version("fftw"),
        module_version("cray-fftw"),
    )
    fiat_library_paths = installed_library_paths(fiat_root)
    field_api_library_paths = installed_library_paths(field_api_root)
    mpi_dependency_version, mpi_version_output = mpi_version(cache)
    lapack_dependency_version = first_value(
        manifest.get("LAPACK_VERSION"),
        lapack_version(fortran_compiler, lapack_library_paths),
    )

    provenance = {
        "system": system_info(),
        "environment": environment_info(),
        "cmake": cmake_info(cache, manifest),
        "compilers": compilers,
        "dependencies": {
            "lapack": dependency(
                "lapack",
                lapack_dependency_vendor,
                lapack_dependency_version,
                lapack_root,
                lapack_library_paths,
            ),
            "fftw": dependency(
                "fftw",
                fftw_dependency_vendor,
                fftw_version,
                fftw_root,
                fftw_library_paths,
                cflags=os.environ.get("ECTRANS_FFTW_CFLAGS"),
            ),
            "mpi": dependency(
                "mpi",
                mpi_vendor(mpi_root, mpi_version_output),
                mpi_dependency_version,
                mpi_root,
                mpi_library_paths,
            ),
            "ecbuild": dependency(
                "ecbuild",
                "ECMWF",
                dependency_version("ecbuild", dependency_refs, cache, manifest),
                ecbuild_root,
                [],
            ),
            "fiat": dependency(
                "fiat",
                "ECMWF",
                dependency_version("fiat", dependency_refs, cache, manifest),
                fiat_root,
                fiat_library_paths,
            ),
            "field_api": dependency(
                "field_api",
                "ECMWF",
                dependency_version("field_api", dependency_refs, cache, manifest),
                field_api_root,
                field_api_library_paths,
            ),
        },
    }

    if args.anonymous:
        anonymize_dependency(provenance["dependencies"]["ecbuild"], ecbuild_root, "${ecbuild_ROOT}")
        anonymize_dependency(provenance["dependencies"]["fiat"], fiat_root, "${fiat_ROOT}")
        anonymize_dependency(provenance["dependencies"]["field_api"], field_api_root, "${field_api_ROOT}")

    secrets = sensitive_env_values()
    provenance = redact(provenance, secrets)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(args.output)


if __name__ == "__main__":
    main()

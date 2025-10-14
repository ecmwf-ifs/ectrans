---
title: Interfacing with Python
---

ecTrans contains a wrapper for accessing a limited subset of functionality from Python, called ectrans4py. Here we present instructions for building ectrans4py and document the available functionality.

## Installing ectrans4py

Installation proceeds at first as with a general installation of ecTrans. Follow [the instructions](installation.html) up to the step "Building ecTrans". Then export `fiat_ROOT` as the full path to the FIAT build directory.

Then, prepare your Python environment. For example, we suggest that you create a new Python venv:

```bash
python3 -m venv venv
. venv/bin/activate
```

Now you can build the ectrans4py wheel:

```bash
pip install build
python -m build --wheel ectrans -o $PWD
pip install ectrans4py*
```

You should now be able to import the ectrans4py module:

```bash
python -c "import ectrans4py"
```

If this command does not work, please [raise an issue](https://github.com/ecmwf-ifs/ectrans/issues) to let us know.

## ectrans4py documentation

Functionality is provided through the `ectransp4y` module. 

### `ectrans_version`

Returns the version string of ecTrans.

_Parameters_:

None

_Returns_:

- **version_string** : _str_  
  The version string.

### `get_legendre_assets`

Fetches the arrays relevant for performing the Legendre transform. Use this function to obtain the Legendre polynomials

_Parameters_:

- 

_Returns_:

### `etrans_inq4py`

### `trans_inq4py`

### `sp2gp_lam4py`

### `gp2sp_lam4py`

### `sp2gp_gauss4py`

### `gp2sp_gauss4py`

### `sp2gp_fft1d4py`

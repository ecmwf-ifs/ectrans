---
title: API
---

# ecTrans API

## General notes

NOTE: ecTrans is a legacy code with an accumulated 30 years of history. Over this time certain
features enabled through optional arguments will have fallen out of use. We are currently reviewing
all options to identify those that can be safely deleted, but this takes time. In the mean time, we
have tagged below all options we deem to be "potentially deprecatable".

### Variable names

ecTrans _in principle_ follows the coding standard and conventions outlined in the [IFS
Documentation - Part VI: Technical and Computational Procedures](https://www.ecmwf.int/en/elibrary/
81372-ifs-documentation-cy48r1-part-vi-technical-and-computational-procedures) section 1.5.
Following these standards, all variable names must begin with a one- or two-character prefix
denoting their scope (module level, dummy argument, local variables, loop index, or parameter) and
type. These are outlined in Table 1.2. Dummy variables have the following prefixes:

- `K` - integer
- `P` - real (single or double precision)
- `LD` - logical
- `CD` - character
- `YD` - derived type

### `KIND` parameters

As with the IFS, integer and real variables in ecTrans always have an explicit `KIND` specification.
These are defined in the [`PARKIND1` module](https://github.com/ecmwf-ifs/fiat/blob/main/src/parkind
/parkind1.F90) which is part of the Fiat library (a dependency of ecTrans). To understand the
subroutines described here, only two must be considered:

- `INTEGER, PARAMETER :: JPIM = SELECTED_INT_KIND(9)` (i.e. 4-byte integer)
- `INTEGER, PARAMETER :: JPRD = SELECTED_REAL_KIND(13,300)` (i.e. 8-byte float)


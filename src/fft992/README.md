
# FFT992 Library

This library contains legacy FFT implementations that were removed from the main ectrans codebase.
It is enabled at compile time by setting the ENABLE_FFT992 cache variable to ON.

## Contents

### FFT992 Algorithm
- **fft992.F90**:    Main FFT992 implementation
- **fft992_cc.F90**: Complex-transform helper with the FFT992 calling interface
- **set99.F90**:     Setup routine for FFT992
- **set99b.F90**:    Alternative setup routine for FFT992

### Bluestein Algorithm
- **bluestein_mod.F90**: Bluestein FFT algorithm implementation

## History

These files were part of the initial ecTrans 1.0.0 release, and removed in the following commits:
- FFT992: commit 1e82391 (Dec 5, 2024) - "removed fft992 sources"
- Bluestein: commit 73ab525 (Jun 4, 2024) - "Remove fft992 and Bluestein FFT functionality"
- GPU Bluestein: commit 441f5d8 (Jan 13, 2022) - "Removal of Bluestein- and FFT992-related files"

## Reason for Removal

These implementations were removed to allow simplification of the codebase and avoiding maintenance of multiple FFT codepaths.

## Purpose of This Library

This library preserves the legacy FFT implementations as a separate component that can be:
- used as a fallback when an FFTW-compatible library is not available
- referenced for algorithmic comparisons
- maintained independently from the main ectrans codebase

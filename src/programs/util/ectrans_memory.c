#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>

#define PREFIX ectrans_memory_
#define CONCAT_(A, B) A##B
#define CONCAT(A, B) CONCAT_(A, B)
#define PREFIXED(symbol) CONCAT(PREFIX, symbol)

#define allocate_var                    PREFIXED(allocate_var)
#define deallocate_var                  PREFIXED(deallocate_var)
#define set_label                       PREFIXED(set_label)
#define unset_label                     PREFIXED(unset_label)
#define set_logging                     PREFIXED(set_logging)
#define set_logging_fortran_output_unit PREFIXED(set_logging_fortran_output_unit)
#define write_to_fortran_unit           PREFIXED(write_to_fortran_unit)
#define set_pinning                     PREFIXED(set_pinning)

// ----------------------------------------------------------------------------------------
// API

// Allocate memory of given bytes
void* allocate_var(size_t bytes);

// Deallocate memory of given bytes
//   The bytes argument is not strictly needed but may be useful for logging
void deallocate_var(void* ptr, size_t bytes);

// Set/Unset a label that can be used for logging during the above calls
void set_label(const char* label);
void unset_label();

// Set logging. accepted values: 1,0 (true, false)
void set_logging(int);

// Set Fortran output unit
void set_logging_fortran_unit(int);

// Set pinning. accepted values: 1,0 (true, false)
// Pinning is off by default
void set_pinning(int);

// ----------------------------------------------------------------------------------------
// Implementation

#if defined(CUDA)
#include <cuda_runtime.h>
#define hicHostRegisterMapped cudaHostRegisterMapped
#define hicHostRegister(ptr, bytes, flag) cudaHostRegister(ptr, bytes, flag);
#define hicHostUnregister(ptr) cudaHostUnregister(ptr);
#define HIC
#elif defined(HIP)
#include <hip/hip_runtime.h>
#define hicHostRegisterMapped hipHostRegisterMapped
#define hicHostRegister(ptr, bytes, flag) hipHostRegister(ptr, bytes, flag);
#define hicHostUnregister(ptr) hipHostUnregister(ptr);
#define HIC
#else
// Mockup version of pinning for host-only
#define hicHostRegisterMapped 0
static int hicHostRegister(void* ptr, size_t bytes, int flag) { return 0; }
static int hicHostUnregister(void* ptr) { return 0; }
#endif

static const char* label_ = NULL;
static int logging_ = 0;
static int logging_fortran_output_unit_ = 0;
static int pinning_ = 0;

void write_to_fortran_unit( int unit, const char* msg );

static void fortran_printf(const char *format, ...) {
    char buffer[256];
    va_list args;
    va_start (args, format);
    vsnprintf (buffer,256,format, args);
    write_to_fortran_unit(logging_fortran_output_unit_, buffer);
    va_end (args);
}

static void default_err_callback(const char* err_msg) {
  fortran_printf(0, "%s\n", err_msg);
}
typedef void (*err_callback_t)(const char* err_msg);
static err_callback_t err_callback = default_err_callback;
// err_callback could be runtime-configured to something that calls e.g. abor1

void set_label(const char* label) {
  label_ = label;
}

void unset_label() {
  label_ = NULL;
}

void set_pinning(int pinning) {
  if (logging_) {
    fortran_printf("Set pinning to %d\n", pinning);
  }
  pinning_ = pinning;
}

void set_logging(int logging) {
  logging_ = logging;
}

void set_logging_fortran_output_unit(int output_unit) {
  logging_fortran_output_unit_ = output_unit;
}

static void* allocate_pinned(size_t bytes) {
  void* ptr = malloc(bytes);
  if (bytes && pinning_) {
    int err = hicHostRegister(ptr, bytes, hicHostRegisterMapped);
    if (err) {
      err_callback("Error pinning memory");
    }
  }
  return ptr;
}

static void deallocate_pinned(void* ptr, size_t bytes) {
  if (ptr && pinning_) {
    int err = hicHostUnregister(ptr);
    if (err) {
        err_callback("Error unpinning memory");
    }
    free(ptr);
  }
}

void* allocate_var(size_t bytes) {
  if (logging_) {
    if (label_) {
      fortran_printf("Allocating variable %s with %zu %s bytes\n", label_, bytes, pinning_ ? "pinned" : "");
    }
    else {
      fortran_printf("Allocating variable with %zu %s bytes\n", bytes, pinning_ ? "pinned" : "");
    }
  }

  return allocate_pinned(bytes);
}

void deallocate_var(void* ptr, size_t bytes) {
  if (logging_) {
    if (label_) {
      fortran_printf("Deallocating variable %s with %zu %s bytes\n", label_, bytes, pinning_ ? "pinned" : "");
    }
    else {
      fortran_printf("Deallocating variable with %zu %s bytes\n", bytes, pinning_ ? "pinned" : "");
    }
  }
  deallocate_pinned(ptr, bytes);
}

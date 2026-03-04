/* (C) Copyright 2026- ECMWF.

This software is licensed under the terms of the Apache Licence Version 2.0
which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
In applying this licence, ECMWF does not waive the privileges and immunities
granted to it by virtue of its status as an intergovernmental organisation
nor does it submit to any jurisdiction. */

// This file takes care of enabling floating point exception trapping for the API tests,
// and defines a signal handler to print a stack trace when an FPE is encountered.
// It was taken and adapted from the github.com/ecmwf/atlas repository,
// where it has been used for a while to enable FPE trapping for the Atlas tests.
//
// Note that on Apple Silicon a SIGFPE may be posing as a SIGILL, so we also check for that case in the signal handler.
//
// This could eventually be moved to a more central location, or be incorporated with DR_HOOK
// if we want to enable FPE trapping for more than just the API tests, but for now this is sufficient.

#include <cfenv>
#include <cstring>
#include <iomanip>
#include <bitset>
#include <iostream>
#include <csignal>
#include <cerrno>

extern "C" void linux_trbk(void);

#ifndef HAVE_FEENABLEEXCEPT
#error "This code requires HAVE_FEENABLEEXCEPT to compile"
#endif
#ifndef HAVE_FEDISABLEEXCEPT
#error "This code requires HAVE_FEDISABLEEXCEPT to compile"
#endif

#if HAVE_FEENABLEEXCEPT
static int ectrans_test_feenableexcept(unsigned int excepts) {
    return ::feenableexcept(excepts);
}
static int ectrans_test_fedisableexcept(unsigned int excepts) {
    return ::fedisableexcept(excepts);
}
#elif defined(__APPLE__)
static int ectrans_test_feenableexcept(unsigned int excepts) {
    static fenv_t fenv;
    unsigned int new_excepts = excepts & FE_ALL_EXCEPT;
    unsigned int old_excepts;   // previous masks

    if (::fegetenv(&fenv)) {
        return -1;
    }

#if defined(__arm64__)
    old_excepts = fenv.__fpsr & FE_ALL_EXCEPT;

    fenv.__fpsr |= new_excepts;
    fenv.__fpcr |= (new_excepts << 8);
#else
    old_excepts = fenv.__control & FE_ALL_EXCEPT;

    fenv.__control &= ~new_excepts;
    fenv.__mxcsr   &= ~(new_excepts << 7);
#endif

    return ::fesetenv(&fenv) ? -1 : old_excepts;
}
static int ectrans_test_fedisableexcept(unsigned int excepts) {
    static fenv_t fenv;
    unsigned int new_excepts = excepts & FE_ALL_EXCEPT;
    unsigned int old_excepts;   // all previous masks

    if (::fegetenv(&fenv)) {
        return -1;
    }

#if defined(__arm64__)
    old_excepts = fenv.__fpsr & FE_ALL_EXCEPT;

    fenv.__fpsr &= ~new_excepts;
    fenv.__fpcr &= ~(new_excepts << 8);
#else
    old_excepts = fenv.__control & FE_ALL_EXCEPT;

    fenv.__control |= new_excepts;
    fenv.__mxcsr   |= (new_excepts << 7);
#endif

    return ::fesetenv(&fenv) ? -1 : old_excepts;
}
#else
static int ectrans_test_feenableexcept(unsigned int excepts) {
    return 0;
}
static int ectrans_test_fedisableexcept(unsigned int excepts) {
    return 0;
}
#endif

[[noreturn]] void ectrans_test_signal_handler(int signum, siginfo_t* si, [[maybe_unused]] void* ucontext) {
    std::string signal_code;
    if (signum == SIGFPE) {
        switch (si->si_code) {
            case FPE_FLTDIV:
                signal_code = " [FE_DIVBYZERO]";
                break;
            case FPE_FLTINV:
                signal_code = " [FE_INVALID]";
                break;
            case FPE_FLTOVF:
                signal_code = " [FE_OVERFLOW]";
                break;
            case FPE_FLTUND:
                signal_code = " [FE_UNDERFLOW]";
                break;
            case FPE_FLTRES:
                signal_code = " [FE_INEXACT]";
                break;
        }
    }
#if defined(__APPLE__) && defined(__arm64__)
    if (signum == SIGILL) {
        // On Apple Silicon a SIGFPE may be posing as a SIGILL
        // See:
        //    https://developer.apple.com/forums/thread/689159?answerId=733736022
        //    https://developer.arm.com/documentation/ddi0595/2020-12/AArch64-Registers/ESR-EL1--Exception-Syndrome-Register--EL1-?lang=en#fieldset_0-24_0_16-1_1
        auto esr = reinterpret_cast<ucontext_t*>(ucontext)->uc_mcontext->__es.__esr;
        auto is_floating_point_exception = [&esr]() {
            constexpr unsigned long fpe_mask = 2952790016; // bits: 10110000000000000000000000000000
            constexpr std::bitset<32> fpe_mask_bits(fpe_mask);
            return((fpe_mask_bits & std::bitset<32>(esr)) == fpe_mask_bits);
        };
        auto test_esr = [&esr](auto pos) -> bool {
            return std::bitset<32>(esr).test(pos);
        };
        if (is_floating_point_exception()) {
            // SIGILL is posing as a SIGFPE
            constexpr size_t IOF = 0; // invalid operation
            constexpr size_t DZF = 1; // divide-by-zero
            constexpr size_t OFF = 2; // overflow
            constexpr size_t UFF = 3; // underflow
            constexpr size_t IXF = 4; // inexact
            constexpr size_t IDF = 7; // denormal
            if (test_esr(IOF)) {
                signal_code = " [FE_INVALID]";
            }
            else if(test_esr(DZF)) {
                signal_code = " [FE_DIVBYZERO]";
            }
            else if(test_esr(OFF)) {
                signal_code = " [FE_OVERFLOW]";
            }
            else if(test_esr(UFF)) {
                signal_code = " [FE_UNDERFLOW]";
            }
            else if(test_esr(IXF)) {
                signal_code = " [FE_INEXACT]";
            }
            else if(test_esr(IDF)) {
                signal_code = " [FE_DENORMAL]";
            }
        }
    }
#endif

    std::ostream& out = std::cerr;
    out << "\n"
        << "=========================================\n"
        << signal_code << " (signal intercepted by ectrans_test_signal_handler [ectrans/tests/trans/api/fpe_trapping.cc])\n"
        << "=========================================\n";
    linux_trbk();
    out << "=========================================\n"
        << std::endl;

    // Restore the default signal handler and re-raise the signal to allow for normal termination and core dump generation.
    std::signal(signum, SIG_DFL);
    std::raise(signum);

    // Just in case we end up here, which normally we shouldn't.
    std::cerr << "Exit\n" << std::endl;
    std::_Exit(EXIT_FAILURE);
}

extern "C" {
void ectrans_test_enable_fpe() {
    char* ECTRANS_TEST_ENABLE_FPE = getenv("ECTRANS_TEST_ENABLE_FPE");
    // Don't enable FPE trapping if the environment variable ECTRANS_TEST_ENABLE_FPE is set to "0",
    // to allow for easier debugging of tests when desired without having to modify the code.
    if (ECTRANS_TEST_ENABLE_FPE != nullptr && std::strcmp(ECTRANS_TEST_ENABLE_FPE, "0") == 0) {
        return;
    }
    struct sigaction sa;
    std::memset(&sa, 0, sizeof(sa));
    sa.sa_sigaction = ectrans_test_signal_handler;
    sa.sa_flags = SA_SIGINFO;

    if (sigaction(SIGFPE, &sa, nullptr) == -1) {
        std::cerr << "Failed to set signal handler for SIGFPE: " << std::strerror(errno) << std::endl;
    }

    if (sigaction(SIGILL, &sa, nullptr) == -1) {
        std::cerr << "Failed to set signal handler for SIGILL: " << std::strerror(errno) << std::endl;
    }

    ectrans_test_feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
}
}
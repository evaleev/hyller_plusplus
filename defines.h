
#ifndef _hyller_defines_h_
#define _hyller_defines_h_

/// Set to 1 to skip Hylleraas calculation
#define SKIP_HYLLERAAS 0
/// Set to 1 to skip Hartree-Fock calculation
#define SKIP_HF 0
/// Set to 1 to skip tests
#define SKIP_TESTS 1
/// Set to 1 to test proper normalization of basis sets. Only useful for testing
#define TEST_NORMALIZATION 1
/// Set to 1 to test correct arguments to fac, etc. Only useful for initial testing
#define CORRECT_RATHER_THAN_SLOW 1

/// Sets the max argument for fac, etc.
#define HYLLER_MAXINDEX 200

#endif

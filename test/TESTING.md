# Testing bwa-mem3

bwa-mem3 tests are organized into three categories — **unit**, **integration**, and **regression** — each with its own build target, CI scope, and runtime cost.

| Category       | Binary / runner                    | Fixtures                                     | CI scope                             |
|----------------|------------------------------------|----------------------------------------------|--------------------------------------|
| **unit**       | `test/bwa_mem3_tests_unit`         | None. All inputs synthetic.                  | Every matrix row                     |
| **integration**| `test/bwa_mem3_tests_integration`  | Small committed FASTAs / FMI under `test/fixtures/` | SSE4.1, AVX2, ARM64 Linux, macOS ARM |
| **regression** | `test/regression/*.sh`             | Downloaded references (phiX, chr22) + bwa + dwgsim | Canonical AVX2 row only      |

Design rationale: see the internal test framework design spec (not committed).

## Running tests locally

All commands run from the repo root unless noted.

### Build everything

```bash
make                                      # builds bwa-mem3 binary
make -C test -j$(nproc 2>/dev/null || sysctl -n hw.ncpu)   # builds framework + both doctest binaries + legacy integration binaries
```

### Run unit tests

```bash
./test/bwa_mem3_tests_unit                # run all unit tests
./test/bwa_mem3_tests_unit --list-test-cases
./test/bwa_mem3_tests_unit --test-case="*kswv*"
./test/bwa_mem3_tests_unit --test-suite="unit/kswv"
./test/bwa_mem3_tests_unit --test-suite-exclude=slow
./test/bwa_mem3_tests_unit --success      # also print passing assertions
./test/bwa_mem3_tests_unit --reporters=junit --out=unit-results.xml
./test/bwa_mem3_tests_unit --help         # full doctest flag list
```

### Run integration tests

```bash
./test/bwa_mem3_tests_integration         # same doctest CLI as unit binary
```

Plus the legacy integration binaries (being migrated to the doctest integration binary one-by-one):

```bash
bash test/run_unit_tests.sh               # runs fmi_test, smem2_test, sa2ref_test, bwt_seed_strategy_test, xeonbsw
```

### Run regression checks

Regression scripts expect upstream CI to have staged their inputs. To run one locally:

```bash
# Example: phiX parity. Needs dwgsim to be installed.
mkdir -p /tmp/ci-test && cd /tmp/ci-test
curl -sL "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/819/615/GCF_000819615.1_ViralProj14015/GCF_000819615.1_ViralProj14015_genomic.fna.gz" \
  | gunzip > phix174.fa
dwgsim -z 42 -N 500 -1 150 -2 150 -r 0.001 -S 2 phix174.fa reads
cd -
BWA_MEM3="$(pwd)/bwa-mem3" CI_TEST_DIR=/tmp/ci-test bash test/regression/phix_parity.sh
```

## Running tests in CI

- **ci.yml** runs unit tests on every matrix row, integration tests on the four widened canonical rows, and regression tests on the canonical AVX2 row only.
- JUnit artifacts are uploaded per row (`unit-results-<name>.xml`, `integration-results-<name>.xml`). Download them from a failed run's Actions page to see fine-grained assertion output.
- **proto-neon-kswv.yml** runs `./test/bwa_mem3_tests_unit --test-suite="unit/kswv"` on proto branches.

## Adding a unit test

Unit tests are doctest `TEST_CASE`s in `test/unit/test_*.cpp`. They must:

1. Use only synthetic inputs (no files on disk).
2. Complete in <100 ms each.
3. Be tagged with `doctest::test_suite("unit/<module>")` where `<module>` is one of `kswv`, `bandedsw`, `ksw`, `fmindex`, `smem`, `bam`, `pair`, `cigar`, `util`. doctest's `test_suite` decorator is overriding (not additive) — chaining `* doctest::test_suite("a") * doctest::test_suite("b")` keeps only the last one — so encode the category and module as a single slash-separated string.

Template:

```cpp
// test/unit/test_<module>.cpp
#include "doctest/doctest.h"
#include "<framework-header>.h"

TEST_CASE("short description of behavior"
          * doctest::test_suite("unit/<module>")) {
    // Arrange
    auto input = /* construct synthetic input */;
    // Act
    auto result = function_under_test(input);
    // Assert
    CHECK(result == expected);
}
```

Use `SUBCASE` for parameterized variants:

```cpp
TEST_CASE("behavior across inputs" * doctest::test_suite("unit/<module>")) {
    SUBCASE("small input")   { /* ... */ }
    SUBCASE("empty input")   { /* ... */ }
    SUBCASE("huge input")    { /* ... */ }
}
```

To filter by category or module, use a glob: `--test-suite="unit/*"` selects everything under `unit/`, `--test-suite="*/kswv"` selects every category's `kswv` tests.

Framework helpers available in `test/framework/`:
- `bwa_tests::build_scoring_matrix`, `default_scoring_matrix` (scoring.h)
- `bwa_tests::gen_{random,exact_match,all_mismatch,homopolymer,sub_cluster,with_n_bases}_pair` (seqpair_gen.h)
- `bwa_tests::BatchBuffers` (seqpair_batch.h)
- `bwa_tests::run_scalar_ksw` (ksw_runner.h)
- `bwa_tests::run_kswv_batch` (kswv_runner.h)
- `bwa_tests::kswr_{score,coords,score2}_eq` (kswr_cmp.h)

The file is picked up automatically by the wildcard rule in `test/Makefile` — no Makefile edits needed when adding a new `test_*.cpp`.

## Adding an integration test

Same as unit tests, but:
- Lives in `test/integration/test_*.cpp` (note the `test_` prefix — plain `.cpp` files in that dir are reserved for the legacy binaries in flight).
- Tagged with `doctest::test_suite("integration/<module>")` (same overriding-decorator caveat as unit).
- May load fixtures from `test/fixtures/`. Fixtures must be committed, ≤100 KB, and come with a documented regeneration procedure.
- Runtime budget per `TEST_CASE`: <10 s.

## Adding a regression test

Write a standalone bash script in `test/regression/<name>.sh`:
- `set -euo pipefail` at the top.
- Document required env-var inputs in the comment block.
- Emit `PASS:` / `FAIL:` lines; exit nonzero on failure.
- Add a row to `test/regression/README.md`.
- Wire it into `.github/workflows/ci.yml` with the appropriate `if: matrix.canonical == true` gate.

Choose regression over integration when:
- The test shells out to bwa-mem3 (or any other binary) as a process.
- The test diffs against a third-party tool's output (bwa, bwa-meth, samtools).
- The test needs fixtures >100 KB or downloaded at CI time.

## Debugging a failing test

```bash
./test/bwa_mem3_tests_unit --test-case="*kswv*" --success  # verbose
./test/bwa_mem3_tests_unit --test-case="*kswv*" --break    # break into debugger on failure
./test/bwa_mem3_tests_unit --test-case="*foo*" --subcase="bar"  # run single SUBCASE

# Re-run kswv tests with the original phase-0/phase-1 diagnostics:
BWA_TESTS_DEBUG_PHASE0=1 BWA_TESTS_DEBUG_PHASE1=1 \
  ./test/bwa_mem3_tests_unit --test-suite="unit/kswv"
```

For a test using `std::mt19937(seed)`, reproduce locally by editing the seed and running the full binary — there is no persistent "last-failing-seed" record.

## Upgrading doctest

1. Pick a new upstream tag from <https://github.com/doctest/doctest/releases>.
2. Fetch the header:
   ```bash
   curl -sSL -o ext/doctest/doctest.h \
     https://raw.githubusercontent.com/doctest/doctest/<new-tag>/doctest/doctest.h
   ```
3. Compute and record the new SHA256 in `ext/doctest/VERSION`:
   ```bash
   shasum -a 256 ext/doctest/doctest.h
   ```
4. Run the full test suite locally on x86 and ARM.
5. Push to a branch and verify CI passes on every matrix row before merging.
6. Commit with message `build(ext): upgrade doctest to <new-tag>`.

## Framework reference

| Header                           | What                                           |
|----------------------------------|------------------------------------------------|
| `framework/scoring.h`            | `ScoringMatrix` type, `build_scoring_matrix`   |
| `framework/seqpair.h`            | `TestPair` struct                              |
| `framework/seqpair_gen.h`        | Deterministic TestPair generators              |
| `framework/seqpair_batch.h`      | `BatchBuffers` — flat-layout packer            |
| `framework/ksw_runner.h`         | Scalar `run_scalar_ksw`, default gap/xtra      |
| `framework/kswv_runner.h`        | Two-pass `run_kswv_batch`                      |
| `framework/kswr_cmp.h`           | Score / coord / score2 comparators             |
| `framework/junit_reporter.h`     | CI matrix-row banner                           |
| `framework/test_main.cpp`        | Shared `main()` (single doctest IMPL TU)       |

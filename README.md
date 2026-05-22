# RDM™ Twin Prime Density Engine

> **Deterministic. Zero Dependencies. Microsecond Speed.**
> The first engine to count and verify twin primes using pure geometric coordinate mapping — no probabilistic methods, no trial division, no black-box libraries.

[![License](https://img.shields.io/badge/license-Sovereign%20Safe%20Zone-blue)](./LICENSE.md)
[![Language](https://img.shields.io/badge/language-Rust-orange)](https://www.rust-lang.org/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19863559.svg)](https://doi.org/10.5281/zenodo.19863559)

---

## What This Is

The **RDM™ Twin Prime Density Engine** counts, predicts, and verifies twin prime distributions using a coordinate-geometry algorithm based on the Lattice Exclusion Identity (LEI). 

Instead of generating primes to divide by, the engine maps integer spaces onto a tightly packed boolean array. It uses an optimized physical traversal to flag composite collisions natively in memory. The un-flagged coordinates remaining after the sweep are guaranteed, mathematically proven twin prime pairs.

This is not a sieve. This is a deterministic structural coordinate mapper.

---

## The Algorithm & Mathematics

This README documents the software engineering, benchmarks, and usage of the engine. 

**For the formal mathematical proofs, Trace-to-Lane Splitting equations, and exact cardinality identities that power this algorithm, please refer to the official manuscripts published on Zenodo. [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19863559.svg)](https://doi.org/10.5281/zenodo.19863559)**

---

## Benchmarks

All benchmarks run on Apple Silicon (release build, `cargo build --release`). No warmup. Cold binary.

| Command | N | Result | Time |
|---|---|---|---|
| `predict` | 1,000,000 | π₂ ≈ 6,917 | **9.8 ms** (includes Euler product over 78,497 primes) |
| `verify` | 100,000 | π₂ = 1,224 (exact) | **223 µs** (LEI scan only) |
| `verify` | 1,000,000 | π₂ = 8,169 (exact) | **2.99 ms** |
| `verify` | 10,000,000 | π₂ = 58,980 (exact) | **33.9 ms** |
| `verify-range` | 1,000 → 100,000 | 1,189 twins, 9,424 primes | **315 µs** |
| `verify-range` | 1,000 → 10,000 | 170 twins in range | **107 µs** |

> The LEI scan at N = 10,000,000 allocated ~1.67 MB of boolean memory (k/6 per lane). No GC. No heap fragmentation.

---

## Installation

Requires [Rust](https://rustup.rs/) (stable, 1.70+).

```bash
git clone https://github.com/trvengine/rdm_twin_density
cd rdm_twin_density
cargo build --release
```

The binary is at `target/release/rdm_twin_density`. Zero runtime dependencies.

---

## Usage

### Show the 2C₂ Twin Prime Constant
```bash
./rdm_twin_density constant
```

### Predict π₂(N)
```bash
./rdm_twin_density predict 1000000
```

### Verify with Exact LEI Count
```bash
./rdm_twin_density verify 1000000
```

### Count Primes and Twins in a Range
```bash
./rdm_twin_density verify-range 1000 10000
```

### Export the Full Coordinate Trace to CSV
```bash
./rdm_twin_density export 100000 trace.csv
```
Output columns: `k, L_minus_prime, L_plus_prime, Is_Twin`

Load this into Python, Excel, or MATLAB to graph the Prime Manifold topology.

---

## Memory Model

The engine allocates exactly two boolean arrays of size $\lfloor N/6 \rfloor + 1$:

```
Memory Used = 2 × ⌊N/6⌋ bytes
```

| N | Memory |
|---|---|
| 1,000,000 | ~333 KB |
| 10,000,000 | ~3.3 MB |
| 100,000,000 | ~33 MB |
| 1,000,000,000 | ~333 MB |

No limit is enforced. The user is responsible for their machine's RAM capacity.

---

## License

This software is released under the **TRV™ Sovereign Safe Zone License**.  
Academic and research use is permitted with attribution.  
Commercial use requires explicit written permission from TRV™ Labs.

See [LICENSE.md](./LICENSE.md) for full terms.

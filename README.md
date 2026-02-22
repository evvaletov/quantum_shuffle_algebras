# Quantum Shuffle Algebras Package

**Mathematica + Julia** | [MIT License](LICENSE.md)

**Author:** Eremey Valetov\
**Creation date:** 09-Sep-2016\
**Email:** eremey@valetov.com

## Description

A package for quantum shuffle algebra calculations, including construction of
bases in terms of Lyndon words (primes), element comparison, unique prime
factorization, quantum shuffle multiplication, Hopf algebra operations
(coproduct, counit, antipode), Lyndon word enumeration, and diagonal braiding.

Available in both **Mathematica** and **Julia**.

## Requirements

- **Mathematica:** version 10.4 or later
- **Julia:** version 1.6 or later

## Installation (Mathematica)

Copy `QuantumShuffleAlgebra.wl` to a directory on Mathematica's search path:

```
$UserBaseDirectory/Applications/
```

Then load:

```mathematica
Get["QuantumShuffleAlgebra`"]
```

## Installation (Julia)

```julia
using Pkg
Pkg.develop(path="julia/")

using QuantumShuffleAlgebras
```

## Quick Start (Mathematica)

```mathematica
Get["QuantumShuffleAlgebra`"]

(* Check whether a word is a Lyndon word *)
QSAIsPrime[{1, 2}]     (* True *)

(* Unique prime factorization *)
QSAUniquePrimeFactorization[{2, 1}]     (* {{2}, {1}} *)

(* Quantum shuffle product *)
QSAShuffleMultiplication[{{1}, {2}}]

(* Enumerate Lyndon words up to length 3 over {1,2} *)
QSALyndonWords[3, {1, 2}]

(* Hopf algebra coproduct *)
QSACoproduct[{1, 2}]
```

## Functions

### Core operations

| Function | Purpose |
|---|---|
| `QSARelation[x, y]` | Compare words (1 = >, 0 = =, -1 = <) |
| `QSAIsPrime[word]` | Test whether a word is a Lyndon word |
| `QSAFirstPrime[word]` | Longest Lyndon prefix |
| `QSAUniquePrimeFactorization[word]` | Factorization into Lyndon words |
| `QSAX[word, qpar]` | Quantum shuffle basis element X_a |
| `QSAShuffleMultiplication[{a, b}]` | Quantum shuffle product v_a . v_b |
| `QSAPrimaryCoefficient[word]` | Coefficient of v_a in X_a |
| `QSAExpressInLyndonWords[word]` | Express v_a in terms of X_c's |

### Lyndon word enumeration

| Function | Purpose |
|---|---|
| `QSALyndonWords[n, alphabet]` | All Lyndon words up to length n |

### Hopf algebra

| Function | Purpose |
|---|---|
| `QSACoproduct[word, q]` | Deconcatenation coproduct |
| `QSACounit[word]` | Counit (1 for empty word, 0 otherwise) |
| `QSAAntipode[word, q]` | Antipode via Hopf axiom |

### Braiding

| Function | Purpose |
|---|---|
| `QSABraid[word1, word2, qMatrix]` | Diagonal braiding |

## Tests

```
wolframscript -file tests/RunTests.wl
```

## Copyright and Citation

© 2016 Eremey Valetov. Released under the [MIT License](LICENSE.md).

This package was developed as part of an internship project at UPMC —
Université Pierre et Marie Curie — Paris 6 (now part of Sorbonne Université),
titled "Bases of Quantum Group Algebras in Terms of Lyndon Words". The project
report is available at <https://hal.science/hal-02448969>
(DOI: [10.48550/arXiv.2001.10435](https://doi.org/10.48550/arXiv.2001.10435)).

If you use this package in your research, please cite:

> Eremey Valetov. Bases of Quantum Group Algebras in Terms of Lyndon Words.
> [Internship report] Université Pierre & Marie Curie — Paris 6; Université
> Paris Diderot — Paris 7. 2016. ⟨hal-02448969v2⟩

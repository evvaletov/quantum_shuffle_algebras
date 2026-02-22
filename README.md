# Quantum Shuffle Algebras Package

**Author:** Eremey Valetov\
**Creation date:** 09-Sep-2016\
**Email:** evv@msu.edu

## Description

A Mathematica package for quantum shuffle algebra calculations, including
construction of bases in terms of Lyndon words (primes), element comparison,
unique prime factorization, and quantum shuffle multiplication.

## Requirements

- Mathematica 10.4 or later

## Installation

Copy `QSA Package.nb` to a directory on Mathematica's search path, for example:

```
$UserBaseDirectory/Applications/
```

You can find this path by evaluating `$UserBaseDirectory` in a Mathematica
notebook. On most systems it resolves to:

- **Linux:** `~/.Mathematica/Applications/`
- **macOS:** `~/Library/Mathematica/Applications/`
- **Windows:** `%APPDATA%\Mathematica\Applications\`

## Quick Start

```mathematica
(* Load the package *)
Get["QSA Package.nb"]

(* Check whether a word is prime (a Lyndon word) *)
QSAIsPrime[{1, 2}]
(* True *)

(* Compute the unique prime factorization of a word *)
QSAUniquePrimeFactorization[{2, 1, 2}]
(* {{1, 2}, {2}} *)
```

## Functions

| Function | Purpose |
|---|---|
| `QSAX[word]` | Represent a word as an element of the algebra |
| `QSARelation[i, j]` | The q-relation between generators i and j |
| `QSAShuffleMultiplication[a, b]` | Quantum shuffle product of two elements |
| `QSAIsPrime[word]` | Test whether a word is a Lyndon word (prime) |
| `QSAFirstPrime[word]` | Extract the first prime factor |
| `QSAUniquePrimeFactorization[word]` | Compute the unique factorization into primes |
| `QSAExpressInLyndonWords[expr]` | Express an element in the Lyndon word basis |

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

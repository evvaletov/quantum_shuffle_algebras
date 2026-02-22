# QuantumShuffleAlgebras.jl

Julia implementation of quantum shuffle algebra operations.

Part of the [quantum_shuffle_algebras](https://github.com/evvaletov/quantum_shuffle_algebras)
project (Mathematica + Julia).

## Installation

```julia
using Pkg
Pkg.develop(path="path/to/quantum_shuffle_algebras/julia")
```

## Usage

```julia
using QuantumShuffleAlgebras

# Lyndon words
is_lyndon(Word([1, 2]))                    # true
lyndon_factorize(Word([2, 1]))             # [Word(2), Word(1)]
lyndon_words(3, [1, 2])                    # 5 Lyndon words

# Shuffle product
quantum_shuffle([1], [2])                  # classical
quantum_shuffle([1], [2]; q=[1 3; 5 1])   # quantum

# Hopf algebra
coproduct(Word([1, 2]))                    # deconcatenation
counit(Word(Int[]))                        # 1
antipode(Word([1]))                        # [([1], -1)]

# Braiding
braid(Word([1]), Word([2]), [1 3; 5 1])    # (Word(2,1), 3)
```

## Tests

```
julia --project -e 'using Pkg; Pkg.test()'
```

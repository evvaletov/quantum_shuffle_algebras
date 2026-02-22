module QuantumShuffleAlgebras

export Word, TensorTerm, TensorExpr
export is_lyndon, first_lyndon_prefix, lyndon_factorize, lyndon_words
export relation, shuffle, quantum_shuffle
export coproduct, counit, antipode
export braid

include("lyndon.jl")
include("shuffle.jl")
include("hopf.jl")
include("braid.jl")

end

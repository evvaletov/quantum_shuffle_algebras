"""
    TensorTerm{T,C}

A single term in a formal linear combination of tensor words:
coefficient * (w1 ⊗ w2 ⊗ ... ⊗ wn), stored as coefficient and a flat word.
"""
struct TensorTerm{T,C}
    word::Vector{T}
    coeff::C
end

"""
    TensorExpr{T,C}

A formal linear combination of tensor products, represented as a list of TensorTerms.
"""
struct TensorExpr{T,C}
    terms::Vector{TensorTerm{T,C}}
end

function Base.show(io::IO, te::TensorExpr)
    if isempty(te.terms)
        print(io, "0")
        return
    end
    for (i, t) in enumerate(te.terms)
        i > 1 && print(io, " + ")
        t.coeff != 1 && print(io, t.coeff, "*")
        print(io, "v_", join(t.word, " ⊗ v_"))
    end
end

"""
    simplify(expr::TensorExpr) -> TensorExpr

Combine like terms and remove zeros.
"""
function simplify(expr::TensorExpr{T,C}) where {T,C}
    d = Dict{Vector{T}, C}()
    for t in expr.terms
        d[t.word] = get(d, t.word, zero(C)) + t.coeff
    end
    terms = [TensorTerm{T,C}(w, c) for (w, c) in d if !iszero(c)]
    sort!(terms, by=t -> t.word)
    TensorExpr{T,C}(terms)
end

"""
    shuffle(w1::Vector{T}, w2::Vector{T}) -> Vector{Tuple{Vector{T}, Int}}

Classical shuffle product of two sequences. Returns list of (permuted_word, coefficient).
All coefficients are 1 for classical shuffle.
"""
function shuffle(w1::Vector{T}, w2::Vector{T}) where T
    isempty(w1) && return [(w2, 1)]
    isempty(w2) && return [(w1, 1)]
    result = Tuple{Vector{T}, Int}[]
    for (w, c) in shuffle(w1[2:end], w2)
        push!(result, (vcat([w1[1]], w), c))
    end
    for (w, c) in shuffle(w1, w2[2:end])
        push!(result, (vcat([w2[1]], w), c))
    end
    return result
end

"""
    quantum_shuffle(w1::Vector{T}, w2::Vector{T}; q=nothing) -> TensorExpr

Quantum shuffle product. If q is nothing, returns the classical shuffle.
If q is a matrix, uses diagonal braiding q[i,j] coefficients.
"""
function quantum_shuffle(w1::Vector{T}, w2::Vector{T};
                         q::Union{Nothing, AbstractMatrix}=nothing) where T
    if q === nothing
        terms = shuffle(w1, w2)
        return simplify(TensorExpr([TensorTerm(w, c) for (w, c) in terms]))
    end

    # Quantum shuffle with braiding matrix
    isempty(w1) && return TensorExpr([TensorTerm(w2, one(eltype(q)))])
    isempty(w2) && return TensorExpr([TensorTerm(w1, one(eltype(q)))])

    C = eltype(q)
    result = TensorTerm{T,C}[]

    # w1[1] * sh(w1[2:end], w2)
    sub1 = quantum_shuffle(w1[2:end], w2; q=q)
    for t in sub1.terms
        push!(result, TensorTerm(vcat([w1[1]], t.word), t.coeff))
    end

    # q-weighted: w2[1] * sh(w1, w2[2:end]) with q-factor
    qfactor = prod(q[w1[i], w2[1]] for i in 1:length(w1))
    sub2 = quantum_shuffle(w1, w2[2:end]; q=q)
    for t in sub2.terms
        push!(result, TensorTerm(vcat([w2[1]], t.word), qfactor * t.coeff))
    end

    return simplify(TensorExpr(result))
end

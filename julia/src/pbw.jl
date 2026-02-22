"""
    standard_factorization(l::Word{T}) -> (Word{T}, Word{T})

Standard factorization of a Lyndon word l into (u, v) where:
- l = uv (concatenation)
- v is the longest proper suffix of l that is a Lyndon word
- Both u and v are Lyndon
"""
function standard_factorization(l::Word{T}) where T
    length(l) <= 1 && error("standard_factorization requires a Lyndon word of length >= 2")
    # Find the longest proper Lyndon suffix
    for k in 2:length(l)
        suffix = l[k:length(l)]
        if is_lyndon(suffix)
            return (l[1:k-1], suffix)
        end
    end
    error("no standard factorization found — is the word Lyndon?")
end

"""
    quantum_bracket(l::Word{T}, Q::AbstractMatrix{C}) -> TensorExpr{T,C}

Compute the quantum root vector [l]_q for a Lyndon word l, defined recursively:
- [eᵢ]_q = eᵢ  (single letter)
- [l]_q = [u]_q ⊗_sh [v]_q - q_factor * [v]_q ⊗_sh [u]_q

where (u,v) is the standard factorization of l and q_factor = ∏ Q[vⱼ, uᵢ]
for all letters vⱼ in v and uᵢ in u.
"""
function quantum_bracket(l::Word{T}, Q::AbstractMatrix{C}) where {T, C}
    if length(l) == 1
        return TensorExpr([TensorTerm(l.letters, one(C))])
    end

    u, v = standard_factorization(l)
    bracket_u = quantum_bracket(u, Q)
    bracket_v = quantum_bracket(v, Q)

    # [u] ⊗_sh [v]
    uv = _shuffle_exprs(bracket_u, bracket_v, Q)
    # [v] ⊗_sh [u]
    vu = _shuffle_exprs(bracket_v, bracket_u, Q)

    # q_factor = ∏ Q[v_j, u_i] for all pairs
    qf = one(C)
    for vj in v
        for ui in u
            qf *= Q[vj, ui]
        end
    end

    # [l] = uv - q_factor * vu
    combined = TensorTerm{T,C}[]
    for t in uv.terms
        push!(combined, TensorTerm(t.word, t.coeff))
    end
    for t in vu.terms
        push!(combined, TensorTerm(t.word, -qf * t.coeff))
    end
    simplify(TensorExpr(combined))
end

"""Shuffle two TensorExprs using quantum shuffle with braiding matrix Q."""
function _shuffle_exprs(a::TensorExpr{T,C}, b::TensorExpr{T,C},
                        Q::AbstractMatrix{C}) where {T, C}
    result = TensorTerm{T,C}[]
    for ta in a.terms
        for tb in b.terms
            sh = quantum_shuffle(ta.word, tb.word; q=Q)
            for ts in sh.terms
                push!(result, TensorTerm(ts.word, ta.coeff * tb.coeff * ts.coeff))
            end
        end
    end
    simplify(TensorExpr(result))
end

"""
    pbw_monomials(degree::Int, alphabet::Vector{T}, Q::AbstractMatrix) -> Vector{TensorExpr}

Enumerate PBW basis monomials up to total degree `degree`.
Each monomial is an ordered product [l₁]^{a₁} ⊗_sh [l₂]^{a₂} ⊗_sh ...
where l₁ > l₂ > ... are Lyndon words (in QSA ordering).

Returns a vector of (TensorExpr, description) pairs where description
is the list of (lyndon_word, exponent) pairs.
"""
function pbw_monomials(degree::Int, alphabet::Vector{T}, Q::AbstractMatrix{C}) where {T, C}
    # Get all Lyndon words up to this degree
    lw = lyndon_words(degree, alphabet)

    # Sort in decreasing order (QSA: shorter-is-greater, so largest first)
    sort!(lw, by=w -> w.letters, lt=(a, b) -> relation(Word(a), Word(b)) == 1)

    # Cache quantum brackets
    brackets = Dict{Word{T}, TensorExpr{T,C}}()
    for l in lw
        brackets[l] = quantum_bracket(l, Q)
    end

    # Generate all multisets of Lyndon words with total degree <= degree
    # Each monomial is a list of (lyndon_word_index, exponent) with total degree = sum of (length * exponent)
    result = Tuple{TensorExpr{T,C}, Vector{Tuple{Word{T}, Int}}}[]

    # The empty monomial (degree 0) is the unit
    push!(result, (TensorExpr([TensorTerm(T[], one(C))]), Tuple{Word{T},Int}[]))

    _generate_pbw_monomials!(result, lw, brackets, Q, degree, 1,
                             TensorExpr([TensorTerm(T[], one(C))]),
                             Tuple{Word{T},Int}[], 0)
    result
end

function _generate_pbw_monomials!(result, lw, brackets, Q::AbstractMatrix{C},
                                  max_degree, start_idx, current_expr,
                                  current_desc, current_degree) where C
    for i in start_idx:length(lw)
        l = lw[i]
        wlen = length(l)
        for exp in 1:div(max_degree - current_degree, wlen)
            new_degree = current_degree + wlen * exp
            new_degree > max_degree && break

            # Compute [l]^exp via iterated shuffle
            power = brackets[l]
            for _ in 2:exp
                power = _shuffle_exprs(power, brackets[l], Q)
            end

            # Multiply current expression by this power
            new_expr = if isempty(current_expr.terms[1].word)
                power
            else
                _shuffle_exprs(current_expr, power, Q)
            end

            new_desc = vcat(current_desc, [(l, exp)])
            push!(result, (new_expr, new_desc))

            # Recurse with higher-indexed (smaller in QSA order) Lyndon words
            _generate_pbw_monomials!(result, lw, brackets, Q, max_degree, i + 1,
                                     new_expr, new_desc, new_degree)
        end
    end
end

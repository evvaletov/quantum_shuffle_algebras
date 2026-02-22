"""
    coproduct(w::Word{T}) -> Vector{Tuple{Word{T}, Word{T}}}

Deconcatenation coproduct: Δ(w) = Σ_{k=0}^{n} w[1:k] ⊗ w[k+1:n].
Returns list of (left_tensor, right_tensor) pairs.
"""
function coproduct(w::Word{T}) where T
    result = Tuple{Word{T}, Word{T}}[]
    for k in 0:length(w)
        left = k == 0 ? Word(T[]) : w[1:k]
        right = k == length(w) ? Word(T[]) : w[k+1:length(w)]
        push!(result, (left, right))
    end
    return result
end

"""
    counit(w::Word) -> Int

Counit: ε(1) = 1, ε(eᵢ) = 0 for any non-empty word.
"""
counit(w::Word) = isempty(w) ? 1 : 0

"""
    antipode(w::Word{T}) -> Vector{Tuple{Vector{T}, Int}}

Antipode S via the recursive Hopf algebra formula:
S(1) = 1, S(w) = -w - Σ_{k=1}^{n-1} S(w[1:k]) ⊗ w[k+1:n].
Returns list of (word, coefficient) pairs.
"""
function antipode(w::Word{T}) where T
    isempty(w) && return [(T[], 1)]
    n = length(w)

    # S(w) = -w - Σ_{k=1}^{n-1} S(w[1:k]) * w[k+1:n]
    result = Dict{Vector{T}, Int}()
    result[w.letters] = get(result, w.letters, 0) - 1

    for k in 1:n-1
        left_terms = antipode(w[1:k])
        right = w[k+1:n].letters
        for (lw, lc) in left_terms
            combined = vcat(lw, right)
            result[combined] = get(result, combined, 0) - lc
        end
    end

    return [(w, c) for (w, c) in result if c != 0]
end

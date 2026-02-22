"""
    braid(w1::Word, w2::Word, Q::AbstractMatrix)

Diagonal braiding: σ(e_{i₁}⊗...⊗e_{iₘ}, e_{j₁}⊗...⊗e_{jₙ}) =
(∏ᵢⱼ Q[iₖ,jₗ]) * (e_{j₁}⊗...⊗e_{jₙ}⊗e_{i₁}⊗...⊗e_{iₘ})

Returns (braided_word, coefficient).
"""
function braid(w1::Word{T}, w2::Word{T}, Q::AbstractMatrix) where T
    coeff = one(eltype(Q))
    for i in w1
        for j in w2
            coeff *= Q[i, j]
        end
    end
    return (vcat(w2, w1), coeff)
end

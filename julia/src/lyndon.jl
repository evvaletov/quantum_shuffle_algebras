"""
    Word{T}

A word in the free monoid — thin wrapper around Vector{T}.
"""
struct Word{T}
    letters::Vector{T}
end

Word(xs::T...) where T = Word(collect(xs))
Base.length(w::Word) = length(w.letters)
Base.getindex(w::Word, i) = w.letters[i]
Base.getindex(w::Word, r::UnitRange) = Word(w.letters[r])
Base.iterate(w::Word, s...) = iterate(w.letters, s...)
Base.eltype(::Type{Word{T}}) where T = T
Base.:(==)(a::Word, b::Word) = a.letters == b.letters
Base.hash(w::Word, h::UInt) = hash(w.letters, h)
Base.show(io::IO, w::Word) = print(io, "Word(", join(w.letters, ","), ")")
Base.isempty(w::Word) = isempty(w.letters)
Base.vcat(a::Word{T}, b::Word{T}) where T = Word(vcat(a.letters, b.letters))

"""
    relation(x::Word, y::Word) -> Int

Compare words lexicographically with shorter-is-greater convention.
Returns 1 (x > y), 0 (x == y), -1 (x < y).
"""
function relation(x::Word, y::Word)
    for i in 1:min(length(x), length(y))
        x[i] < y[i] && return -1
        x[i] > y[i] && return 1
    end
    length(x) < length(y) && return 1
    length(x) > length(y) && return -1
    return 0
end

"""
    is_lyndon(w::Word) -> Bool

Test whether a word is a Lyndon word (prime) in the QSA ordering
(shorter-is-greater). A word w is prime if every proper suffix is
strictly less than w, i.e., relation(suffix, w) == -1.
"""
function is_lyndon(w::Word)
    length(w) <= 1 && return true
    for i in 2:length(w)
        relation(Word(w.letters[i:end]), w) > -1 && return false
    end
    return true
end

"""
    first_lyndon_prefix(w::Word) -> Word

Find the longest prefix of w that is a Lyndon word.
"""
function first_lyndon_prefix(w::Word{T}) where T
    best = Word(T[])
    for i in 1:length(w)
        candidate = w[1:i]
        is_lyndon(candidate) || break
        best = candidate
    end
    return best
end

"""
    lyndon_factorize(w::Word) -> Vector{Word}

Compute the unique factorization of w into non-increasing Lyndon words.
"""
function lyndon_factorize(w::Word{T}) where T
    factors = Word{T}[]
    remaining = w
    while !isempty(remaining)
        p = first_lyndon_prefix(remaining)
        push!(factors, p)
        if length(p) < length(remaining)
            remaining = remaining[length(p)+1:length(remaining)]
        else
            remaining = Word(T[])
        end
    end
    return factors
end

"""
    lyndon_words(n::Int, alphabet) -> Vector{Word}

Enumerate all Lyndon words up to length n over the given alphabet.
"""
function lyndon_words(n::Int, alphabet::Vector{T}) where T
    sorted = sort(alphabet)
    words = Word{T}[]
    for a in sorted
        push!(words, Word([a]))
    end
    for len in 2:n
        # Generate all words of this length and filter for Lyndon
        prev = Iterators.product(ntuple(_ -> sorted, len)...)
        for combo in prev
            w = Word(collect(combo))
            is_lyndon(w) && push!(words, w)
        end
    end
    return words
end

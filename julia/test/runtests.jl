using Test
using QuantumShuffleAlgebras

@testset "QuantumShuffleAlgebras" begin

    @testset "Word basics" begin
        w = Word([1, 2, 3])
        @test length(w) == 3
        @test w[1] == 1
        @test w[2:3] == Word([2, 3])
        @test !isempty(w)
        @test isempty(Word(Int[]))
    end

    @testset "Relation (lexicographic ordering)" begin
        @test relation(Word([1]), Word([2])) == -1
        @test relation(Word([2]), Word([1])) == 1
        @test relation(Word([1]), Word([1])) == 0
        # Shorter is greater (QSA convention)
        @test relation(Word([1]), Word([1, 2])) == 1
        @test relation(Word([1, 2]), Word([1])) == -1
    end

    @testset "Lyndon word tests (QSA ordering)" begin
        # In QSA ordering (shorter-is-greater), Lyndon = every suffix < word
        # {2,1} is Lyndon: suffix {1} < {2,1} (1 < 2 at pos 1)
        @test is_lyndon(Word([2, 1])) == true
        # {1,2} is NOT Lyndon: suffix {2} > {1,2} (2 > 1 at pos 1)
        @test is_lyndon(Word([1, 2])) == false
        @test is_lyndon(Word([1])) == true
        @test is_lyndon(Word([1, 1])) == false
        # {3,1,2} is Lyndon: suffix {1,2} < {3,1,2} (1<3), suffix {2} < {3,1,2} (2<3)
        @test is_lyndon(Word([3, 1, 2])) == true
        @test is_lyndon(Word([2, 1])) == true
    end

    @testset "First Lyndon prefix" begin
        @test first_lyndon_prefix(Word([2, 1])) == Word([2, 1])
        @test first_lyndon_prefix(Word([1, 2])) == Word([1])
    end

    @testset "Lyndon factorization" begin
        @test lyndon_factorize(Word([2, 1])) == [Word([2, 1])]
        @test lyndon_factorize(Word([1, 2])) == [Word([1]), Word([2])]
        # {18, 19, 4, 8, 5, 7} → {{18}, {19, 4, 8, 5, 7}}
        factors = lyndon_factorize(Word([18, 19, 4, 8, 5, 7]))
        @test length(factors) == 2
        @test factors[1] == Word([18])
        @test factors[2] == Word([19, 4, 8, 5, 7])
    end

    @testset "Lyndon word enumeration" begin
        # QSA Lyndon words of length <= 2 over {1,2}: {1}, {2}, {2,1}
        lw = lyndon_words(2, [1, 2])
        @test length(lw) == 3
        @test Word([1]) in lw
        @test Word([2]) in lw
        @test Word([2, 1]) in lw

        # Length <= 3 over {1,2}: {1}, {2}, {2,1}, {2,1,1}, {2,2,1} = 5
        lw3 = lyndon_words(3, [1, 2])
        @test length(lw3) == 5

        # All should be Lyndon
        @test all(is_lyndon, lyndon_words(4, [1, 2]))
    end

    @testset "Classical shuffle" begin
        result = shuffle([1], [2])
        words = [r[1] for r in result]
        @test [1, 2] in words
        @test [2, 1] in words
        @test length(result) == 2

        @test shuffle(Int[], [1, 2]) == [([1, 2], 1)]
    end

    @testset "Quantum shuffle" begin
        expr = quantum_shuffle([1], [2])
        @test length(expr.terms) == 2

        Q = [1 3; 5 1]
        expr_q = quantum_shuffle([1], [2]; q=Q)
        @test length(expr_q.terms) == 2
    end

    @testset "Coproduct" begin
        cp = coproduct(Word([1]))
        @test length(cp) == 2
        @test (Word(Int[]), Word([1])) in cp
        @test (Word([1]), Word(Int[])) in cp

        cp2 = coproduct(Word([1, 2]))
        @test length(cp2) == 3
    end

    @testset "Counit" begin
        @test counit(Word(Int[])) == 1
        @test counit(Word([1])) == 0
        @test counit(Word([1, 2])) == 0
    end

    @testset "Antipode" begin
        s1 = antipode(Word([1]))
        @test length(s1) == 1
        @test s1[1] == ([1], -1)

        s0 = antipode(Word(Int[]))
        @test s0 == [(Int[], 1)]
    end

    @testset "Braiding" begin
        Q = [1 3; 5 1]
        w, c = braid(Word([1]), Word([2]), Q)
        @test w == vcat(Word([2]), Word([1]))
        @test c == 3  # Q[1,2]

        w2, c2 = braid(Word([1, 2]), Word([2]), Q)
        @test c2 == Q[1, 2] * Q[2, 2]  # 3 * 1 = 3
    end

    @testset "Standard factorization" begin
        # {2,1} is Lyndon with standard factorization ({2}, {1})
        # {1} is the longest proper Lyndon suffix
        u, v = standard_factorization(Word([2, 1]))
        @test u == Word([2])
        @test v == Word([1])

        # {3,1,2} = ({3}, {1,2})? Let's check: suffix {2} is Lyndon (length 1),
        # suffix {1,2} — is it Lyndon? is_lyndon({1,2}) = false (suffix {2} > {1,2})
        # So longest proper Lyndon suffix is {2}
        u2, v2 = standard_factorization(Word([3, 1, 2]))
        @test is_lyndon(v2)
        @test is_lyndon(u2)
        @test vcat(u2, v2) == Word([3, 1, 2])
    end

    @testset "Quantum bracket" begin
        Q = [1 2; 3 1]

        # Single letter: bracket is just the letter
        b1 = quantum_bracket(Word([1]), Q)
        @test length(b1.terms) == 1
        @test b1.terms[1].word == [1]
        @test b1.terms[1].coeff == 1

        # {2,1}: standard_factorization = ({2}, {1})
        # q_factor = Q[1,2] = 2
        # [2,1]_Q = sh_Q([2],[1]) - 2*sh_Q([1],[2])
        # sh_Q([2],[1]) = {[2,1]=>1, [1,2]=>Q[2,1]} = {[2,1]=>1, [1,2]=>3}
        # sh_Q([1],[2]) = {[1,2]=>1, [2,1]=>Q[1,2]} = {[1,2]=>1, [2,1]=>2}
        # Result = {[2,1]=>1-4, [1,2]=>3-2} = {[2,1]=>-3, [1,2]=>1}
        b21 = quantum_bracket(Word([2, 1]), Q)
        @test !isempty(b21.terms)
        @test length(b21.terms) == 2
        d = Dict(t.word => t.coeff for t in b21.terms)
        @test d[[2, 1]] == -3
        @test d[[1, 2]] == 1

        # With Q = all ones, commutator vanishes (classical case)
        Q_trivial = [1 1; 1 1]
        b21_trivial = quantum_bracket(Word([2, 1]), Q_trivial)
        @test isempty(b21_trivial.terms)
    end

    @testset "PBW monomials" begin
        Q = [1 1; 1 1]

        # Degree 1 over {1,2}: Lyndon words are {1} and {2}
        # PBW monomials: unit, [1], [2] (degree 0 and 1)
        monos = pbw_monomials(1, [1, 2], Q)
        # Should have: (), (1), (2)
        @test length(monos) >= 3

        # Degree 2: adds [1]^2, [2]^2, [1]*[2], and [2,1]
        monos2 = pbw_monomials(2, [1, 2], Q)
        @test length(monos2) > length(monos)
    end

end

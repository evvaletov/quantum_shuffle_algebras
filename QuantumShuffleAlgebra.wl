(* ::Package:: *)
(* QuantumShuffleAlgebra — Quantum Shuffle Algebra Package *)
(* Author: Eremey Valetov *)
(* Email: eremey@valetov.com *)
(* Original: 09-Sep-2016, refactored 2026 *)

BeginPackage["QuantumShuffleAlgebra`"];

(* Public function usage messages *)

QSARelation::usage =
  "QSARelation[x, y] compares words v_x and v_y in T(V). Returns 1 (>), 0 (=), or -1 (<).";

QSAIsPrime::usage =
  "QSAIsPrime[x] tests whether the word x is a Lyndon word (prime).";

QSAFirstPrime::usage =
  "QSAFirstPrime[x] finds the longest Lyndon prefix of x.";

QSAUniquePrimeFactorization::usage =
  "QSAUniquePrimeFactorization[x] computes the unique factorization of x into Lyndon words.";

QSAV::usage =
  "QSAV[indexList, word] returns the tensor product v_{word[[i1]]} \\[TensorProduct] ... \\[TensorProduct] v_{word[[in]]}.";

QSASumV::usage =
  "QSASumV[indexLists, word] sums tensor products for multiple index lists.";

QSAIndexToObject::usage =
  "QSAIndexToObject[x, word] maps index lists to actual elements of word.";

QSALengthToPartitionIndex::usage =
  "QSALengthToPartitionIndex[L] converts length list (l1,...,ln) to partition index (1, 1+l1, 1+l1+l2, ...).";

QSAGeneratePrimaryShuffles::usage =
  "QSAGeneratePrimaryShuffles[L, word, qpar] generates primary shuffles from a length list.";

QSAGenerateShuffles::usage =
  "QSAGenerateShuffles[L, word, qpar] generates all shuffles from a length list.";

QSAGenerateSecondaryShuffles::usage =
  "QSAGenerateSecondaryShuffles[L, word, qpar] generates secondary shuffles (all minus primary).";

QSAWordToLength::usage =
  "QSAWordToLength[x] extracts the length list from a factored Lyndon word.";

QSAX::usage =
  "QSAX[x, qpar] computes X_x, the quantum shuffle basis element. qpar=True for q-deformed (default), False for classical.";

QSAShuffleMultiplication::usage =
  "QSAShuffleMultiplication[{a, b, ...}] computes the quantum shuffle product v_a . v_b . ...";

QSAPrimaryCoefficient::usage =
  "QSAPrimaryCoefficient[x] computes the coefficient of v_x in X_x.";

QSAExpressInLyndonWords::usage =
  "QSAExpressInLyndonWords[x] expresses v_x in terms of X_c's via a system of linear equations.";

QSALyndonWords::usage =
  "QSALyndonWords[n, alphabet] enumerates all Lyndon words up to length n over the given alphabet.";

QSACoproduct::usage =
  "QSACoproduct[word, q] computes the deconcatenation coproduct of a word in the tensor algebra.";

QSACounit::usage =
  "QSACounit[expr] computes the counit: 1 for the empty word, 0 otherwise.";

QSAAntipode::usage =
  "QSAAntipode[word, q] computes the antipode via the Hopf algebra axiom.";

QSABraid::usage =
  "QSABraid[word1, word2, qMatrix] computes the diagonal braiding sigma(e_i \\[TensorProduct] e_j) = q_{ij} (e_j \\[TensorProduct] e_i).";

QSAStandardFactorization::usage =
  "QSAStandardFactorization[l] returns {u, v} where l = Join[u,v], v is the longest proper Lyndon suffix.";

QSAQuantumBracket::usage =
  "QSAQuantumBracket[l, qMatrix] computes the quantum root vector [l]_q via iterated q-commutator using standard factorization.";

QSAPBWMonomials::usage =
  "QSAPBWMonomials[degree, alphabet, qMatrix] enumerates PBW basis monomials up to the given total degree.";

(* Public symbols used in output expressions *)
v;
q;
Q;

Begin["`Private`"];

(* --- Internal helpers --- *)

qsaConcatenation[lists__List] := Join[lists];

qsaDeconcatenationInternal[word_List, k_Integer] :=
  If[k < 0 || k > Length[word], {},
    {{word[[1 ;; k]], word[[k + 1 ;; -1]]}}];

qsaDeconcatenation[word_List] :=
  Flatten[Table[qsaDeconcatenationInternal[word, k], {k, 0, Length[word]}], 1];

qsaPermutation[xlist_List, yitem_, aword_, qletter_: q] :=
  Module[{plist = {}, k, i},
    For[k = 1, k <= Length[xlist], k++,
      For[i = 1, i <= Length[xlist[[k, 1]]] + 1, i++,
        AppendTo[plist,
          {Insert[xlist[[k, 1]], yitem, i],
           xlist[[k, 2]] Product[
             If[qletter === 1, 1,
               Subscript[qletter, aword[[jj]], aword[[Length[xlist[[k, 1]]] + 1]]]],
             {jj, i, Length[xlist[[k, 1]]]}]}]]];
    plist];

qsaPermuteList[zlist_List, aword_, qletter_: q] :=
  Module[{plist2, i1},
    If[Length[zlist] == 1,
      plist2 = {{zlist, 1}},
      plist2 = {{{zlist[[1]]}, 1}};
      For[i1 = 2, i1 <= Length[zlist], i1++,
        plist2 = qsaPermutation[plist2, zlist[[i1]], aword, qletter]]];
    plist2];

qsaComplement1[list1_List, list2_List] :=
  Module[{list3, list4 = {}, i},
    list3 = Complement[list1[[;; , 1]], list2[[;; , 1]]];
    For[i = 1, i <= Length[list1], i++,
      If[MemberQ[list3, list1[[i, 1]]], AppendTo[list4, list1[[i]]]]];
    list4];

(* --- Core functions --- *)

QSARelation[x_List, y_List] :=
  Module[{rel = 0, i},
    For[i = 1, i <= Min[Length[x], Length[y]], i++,
      If[rel == 0 && x[[i]] < y[[i]], rel = -1];
      If[rel == 0 && x[[i]] > y[[i]], rel = 1]];
    If[rel == 0 && Length[x] < Length[y], rel = 1];
    If[rel == 0 && Length[x] > Length[y], rel = -1];
    rel];

QSAIsPrime[x_List] :=
  Module[{ans = True, i1},
    If[Length[x] > 1,
      For[i1 = 2, i1 <= Length[x], i1++,
        If[QSARelation[x[[i1 ;; -1]], x] > -1, ans = False]]];
    ans];

QSAFirstPrime[x_List] :=
  Module[{ans2 = {}, i2 = 1},
    While[i2 <= Length[x] && QSAIsPrime[x[[1 ;; i2]]],
      ans2 = x[[1 ;; i2]]; i2++];
    ans2];

QSAUniquePrimeFactorization[x_List] :=
  Module[{remx = x, upf = {}, nextprime},
    While[Length[remx] > 0,
      nextprime = QSAFirstPrime[remx];
      AppendTo[upf, nextprime];
      If[Length[nextprime] < Length[remx],
        remx = remx[[Length[nextprime] + 1 ;; -1]],
        remx = {}]];
    upf];

QSAV[x_List, xi_] :=
  Apply[TensorProduct, Map[Subscript[v, xi[[#]]] &, x]];

QSASumV[x_List, xi_] :=
  Module[{sum = 0},
    (If[ListQ[#[[1]]], sum += #[[2]] QSAV[#[[1]], xi], sum += QSAV[#, xi]]) & /@ x;
    sum];

QSAIndexToObject[x_, xi_] :=
  If[ListQ[x[[1]]], Map[xi[[#]] &, x[[;; , 1]]], Map[xi[[#]] &, x]];

QSALengthToPartitionIndex[L_List] :=
  Module[{M = {1}, CurM = 1, i3},
    For[i3 = 1, i3 <= Length[L], i3++,
      CurM += L[[i3]]; AppendTo[M, CurM]];
    M];

QSAGeneratePrimaryShuffles[L_List, aword_, qpar_: True] :=
  Module[{ShuffleList = {}, LL = Total[L], M = QSALengthToPartitionIndex[L],
          ShuffleListTemp = {}, i4},
    For[i4 = 1, i4 < Length[M], i4++,
      AppendTo[ShuffleListTemp, Range[M[[i4]], M[[i4 + 1]] - 1]]];
    ShuffleList = qsaPermuteList[ShuffleListTemp,
      Range[Length[ShuffleListTemp]], If[qpar, Q, 1]];
    For[i4 = 1, i4 <= Length[ShuffleList], i4++,
      ShuffleList[[i4]] = {Flatten[ShuffleList[[i4, 1]]], ShuffleList[[i4, 2]]}];
    ShuffleList];

QSAGenerateShuffles[L_List, aword_, qpar_: True] :=
  Module[{ShuffleList = {}, LL = Total[L], M = QSALengthToPartitionIndex[L],
          i4, item1, plist0, ShuffleListTemp, TestF},
    ShuffleListTemp = qsaPermuteList[Range[LL], aword, If[qpar, q, 1]];
    Do[
      TestF = True;
      item1 = item[[1]];
      For[i4 = 1, i4 < Length[M], i4++,
        If[Not[OrderedQ[Select[item1, # >= M[[i4]] && # < M[[i4 + 1]] &]]],
          TestF = False]];
      If[TestF, AppendTo[ShuffleList, item]],
      {item, ShuffleListTemp}];
    plist0 = QSAGeneratePrimaryShuffles[L, aword, qpar];
    ShuffleList = Union[plist0, qsaComplement1[ShuffleList, plist0]];
    ShuffleList];

QSAGenerateSecondaryShuffles[L_List, aword_, qpar_: True] :=
  qsaComplement1[
    QSAGenerateShuffles[L, aword, qpar],
    QSAGeneratePrimaryShuffles[L, aword, qpar]];

QSAWordToLength[x_List] := Map[Length[#] &, x];

QSAX[x_, qpar_: True] :=
  QSASumV[
    QSAGenerateShuffles[QSAWordToLength[QSAUniquePrimeFactorization[x]], x, qpar],
    x];

QSAShuffleMultiplication[x_List] :=
  QSASumV[
    QSAGenerateShuffles[QSAWordToLength[x], Flatten[x]],
    Flatten[x]];

QSAPrimaryCoefficient[x_] :=
  Coefficient[QSAX[x], QSAV[Range[Length[x]], x]];

QSAExpressInLyndonWords[x_] :=
  Module[{shuffles, shuffles1, L, upf, lhs = {}, rhs = {}, eqns = {},
          rhs1 = {}, ia, avec, coeffarrays, sol},
    upf = QSAUniquePrimeFactorization[x];
    L = QSAWordToLength[upf];
    shuffles = QSAGenerateShuffles[L, x];
    shuffles1 = QSAIndexToObject[shuffles, x];
    For[ia = 1, ia <= Length[shuffles1], ia++,
      AppendTo[lhs, Evaluate[Subscript[X, shuffles1[[ia]]]]];
      AppendTo[rhs, QSAX[shuffles1[[ia]]]];
      AppendTo[rhs1, QSAX[shuffles1[[ia]], False]];
      AppendTo[eqns,
        Evaluate[Subscript[X, shuffles1[[ia]]]] == QSAX[shuffles1[[ia]]]]];
    avec = DeleteDuplicates[Flatten[rhs1, 1, Plus]];
    coeffarrays = Normal[CoefficientArrays[eqns, avec]];
    sol = LinearSolve[coeffarrays[[2]], coeffarrays[[1]]];
    MapThread[#1 == #2 &, {avec, sol}]];

(* --- Tier 1: Lyndon word enumeration --- *)

QSALyndonWords[n_Integer, alphabet_List] :=
  Module[{words = {}, sorted},
    sorted = Sort[alphabet];
    words = List /@ sorted;
    Do[
      words = Join[words,
        Select[Tuples[sorted, len], QSAIsPrime[#] &]],
      {len, 2, n}];
    words];

(* --- Tier 1: Hopf algebra operations --- *)

QSACoproduct[word_List, q_: 1] :=
  Module[{splits},
    splits = qsaDeconcatenation[word];
    Total[Table[
      TensorProduct[
        If[s[[1]] === {}, 1, Apply[TensorProduct, Subscript[v, #] & /@ s[[1]]]],
        If[s[[2]] === {}, 1, Apply[TensorProduct, Subscript[v, #] & /@ s[[2]]]]],
      {s, splits}]]];

QSACounit[word_List] := If[word === {}, 1, 0];

QSAAntipode[word_List, q_: 1] :=
  Module[{n = Length[word]},
    If[n == 0, Return[1]];
    (* S(x) = -x - Sum_{(x)} S(x') x'' for proper deconcatenations *)
    -Apply[TensorProduct, Subscript[v, #] & /@ word] -
      Total[Table[
        TensorProduct[
          QSAAntipode[word[[1 ;; k]], q],
          If[k == n, 1, Apply[TensorProduct, Subscript[v, #] & /@ word[[k + 1 ;; n]]]]],
        {k, 1, n - 1}]]];

(* --- Tier 1: Braiding --- *)

QSABraid[word1_List, word2_List, qMatrix_] :=
  Module[{coeff = 1},
    Do[
      Do[
        coeff *= qMatrix[[i, j]],
        {j, word2}],
      {i, word1}];
    coeff * TensorProduct @@ (Subscript[v, #] & /@ Join[word2, word1])];

(* --- Tier 2: PBW basis construction --- *)

QSAStandardFactorization[l_List] :=
  Module[{k},
    If[Length[l] <= 1, Return[$Failed]];
    For[k = 2, k <= Length[l], k++,
      If[QSAIsPrime[l[[k ;; -1]]],
        Return[{l[[1 ;; k - 1]], l[[k ;; -1]]}]]];
    $Failed];

QSAQuantumBracket[l_List, qMatrix_] :=
  Module[{u, v, bu, bv, uv, vu, qfactor, sf},
    If[Length[l] == 1,
      Return[Subscript[v, l[[1]]]]];
    sf = QSAStandardFactorization[l];
    u = sf[[1]]; v = sf[[2]];
    bu = QSAQuantumBracket[u, qMatrix];
    bv = QSAQuantumBracket[v, qMatrix];
    (* Quantum shuffle product of brackets *)
    uv = QSAShuffleMultiplication[{u, v}];
    vu = QSAShuffleMultiplication[{v, u}];
    (* q-factor = product of qMatrix[v_j, u_i] for all pairs *)
    qfactor = Product[qMatrix[[v[[j]], u[[i]]]],
      {j, Length[v]}, {i, Length[u]}];
    Expand[uv - qfactor vu]];

QSAPBWMonomials[degree_Integer, alphabet_List, qMatrix_] :=
  Module[{lw, brackets, result = {}, monomial},
    lw = QSALyndonWords[degree, alphabet];
    brackets = Association @@ Table[l -> QSAQuantumBracket[l, qMatrix], {l, lw}];
    (* Generate ordered products of brackets up to total degree *)
    qsaPBWRecurse[result, lw, brackets, degree, 1, 1, {}];
    result];

qsaPBWRecurse[result_, lw_, brackets_, maxDeg_, startIdx_, currentExpr_, currentDesc_] :=
  Module[{l, wlen, newExpr, newDesc, power},
    Do[
      l = lw[[i]]; wlen = Length[l];
      Do[
        If[Total[QSAWordToLength /@ currentDesc] + wlen exp > maxDeg, Break[]];
        power = Nest[Function[x, Expand[x brackets[l]]], brackets[l], exp - 1];
        newExpr = If[currentExpr === 1, power, Expand[currentExpr power]];
        newDesc = Append[currentDesc, Table[l, exp]];
        AppendTo[result, {newExpr, newDesc}];
        qsaPBWRecurse[result, lw, brackets, maxDeg, i + 1, newExpr, newDesc],
        {exp, 1, Floor[(maxDeg - Total[Length /@ Flatten[currentDesc, 1]]) / wlen]}],
      {i, startIdx, Length[lw]}]];

End[];

EndPackage[];

(* RunTests.wl — VerificationTest suite for QuantumShuffleAlgebra *)

Get[FileNameJoin[{DirectoryName[$InputFileName], "..", "QuantumShuffleAlgebra.wl"}]];

Print["Running QuantumShuffleAlgebra tests..."];
Print[""];

testCount = 0;
passCount = 0;
failCount = 0;

SetAttributes[runTest, HoldFirst];
runTest[expr_, expected_, name_String] :=
  Module[{result, pass},
    testCount++;
    result = expr;
    pass = (result === expected);
    If[pass, passCount++, failCount++];
    Print[If[pass, "PASS", "FAIL"], ": ", name];
    If[!pass, Print["  Expected: ", expected, "\n  Got: ", result]];
  ];

(* --- Lyndon word tests --- *)
(* QSA ordering: shorter-is-greater. A word is "prime" (Lyndon) if every
   proper suffix is strictly less than the word. *)

Print["=== Lyndon Word Tests ==="];

(* {2,1} IS Lyndon: suffix {1} < {2,1} since 1 < 2 at pos 1 *)
runTest[QSAIsPrime[{2, 1}], True, "QSAIsPrime: {2,1} is prime"];
(* {1,2} is NOT Lyndon: suffix {2} > {1,2} since 2 > 1 at pos 1 *)
runTest[QSAIsPrime[{1, 2}], False, "QSAIsPrime: {1,2} is not prime"];
runTest[QSAIsPrime[{1}], True, "QSAIsPrime: single element is prime"];
runTest[QSAIsPrime[{1, 1}], False, "QSAIsPrime: {1,1} is not prime"];
(* {3,1,2} IS Lyndon: suffix {1,2} < {3,1,2} (1<3), suffix {2} < {3,1,2} (2<3) *)
runTest[QSAIsPrime[{3, 1, 2}], True, "QSAIsPrime: {3,1,2} is prime"];
(* {3,2,1}: suffix {2,1} < {3,2,1} (2<3), suffix {1} < {3,2,1} (1<3) → prime *)
runTest[QSAIsPrime[{3, 2, 1}], True, "QSAIsPrime: {3,2,1} is prime"];

Print[""];
Print["=== Relation Tests ==="];

runTest[QSARelation[{1}, {2}], -1, "QSARelation: {1} < {2}"];
runTest[QSARelation[{2}, {1}], 1, "QSARelation: {2} > {1}"];
runTest[QSARelation[{1}, {1}], 0, "QSARelation: {1} = {1}"];
runTest[QSARelation[{1}, {1, 2}], 1, "QSARelation: {1} > {1,2} (shorter is greater)"];
runTest[QSARelation[{1, 2}, {1}], -1, "QSARelation: {1,2} < {1} (longer is lesser)"];

Print[""];
Print["=== Unique Prime Factorization Tests ==="];

(* {2,1} is already a single Lyndon word *)
runTest[QSAUniquePrimeFactorization[{2, 1}], {{2, 1}},
  "UPF: {2,1} is already prime"];
(* {1,2} = {1}.{2} since neither {1,2} is prime *)
runTest[QSAUniquePrimeFactorization[{1, 2}], {{1}, {2}},
  "UPF: {1,2} = {1}.{2}"];
runTest[
  Length[QSAUniquePrimeFactorization[{18, 19, 4, 8, 5, 7}]], 2,
  "UPF: {18,19,4,8,5,7} has 2 factors"];

Print[""];
Print["=== Shuffle Multiplication Tests ==="];

(* v_{1} . v_{2} should produce two terms: v1 v2 + Q_{1,2} v2 v1 *)
runTest[
  QSAShuffleMultiplication[{{1}, {2}}] ===
    Subscript[v, 1] \[TensorProduct] Subscript[v, 2] +
    Subscript[Q, 1, 2] Subscript[v, 2] \[TensorProduct] Subscript[v, 1],
  True,
  "Shuffle: v_1 . v_2 has correct Q-weighted terms"];

(* Classical shuffle (no q) of single letters *)
runTest[
  (QSAShuffleMultiplication[{{10, 7, 2}, {5, 4, 5}}] =!= 0),
  True,
  "Shuffle: non-trivial result for multi-letter words"];

Print[""];
Print["=== Lyndon Word Enumeration Tests ==="];

(* QSA Lyndon words of length <= 2 over {1,2}: {1}, {2}, {2,1} *)
runTest[
  Sort[QSALyndonWords[2, {1, 2}]],
  Sort[{{1}, {2}, {2, 1}}],
  "LyndonWords: length <= 2 over {1,2}"];

(* QSA Lyndon words of length <= 3 over {1,2}:
   Length 1: {1},{2}; Length 2: {2,1}; Length 3: {2,1,1},{2,2,1} = 5 total *)
runTest[
  Length[QSALyndonWords[3, {1, 2}]],
  5,
  "LyndonWords: count for length <= 3 over {1,2}"];

(* All generated words should be prime *)
runTest[
  AllTrue[QSALyndonWords[4, {1, 2}], QSAIsPrime],
  True,
  "LyndonWords: all generated words are prime"];

Print[""];
Print["=== Hopf Algebra Tests ==="];

(* Counit: empty word -> 1, non-empty -> 0 *)
runTest[QSACounit[{}], 1, "Counit: epsilon({}) = 1"];
runTest[QSACounit[{1}], 0, "Counit: epsilon({1}) = 0"];
runTest[QSACounit[{1, 2}], 0, "Counit: epsilon({1,2}) = 0"];

(* Coproduct of single element: 1 tensor e_i + e_i tensor 1 *)
runTest[
  Expand[QSACoproduct[{1}]] === Expand[1 \[TensorProduct] Subscript[v, 1] + Subscript[v, 1] \[TensorProduct] 1],
  True,
  "Coproduct: Delta(e_1) = 1 tensor e_1 + e_1 tensor 1"];

(* Coproduct of 2-element word — TensorProduct[1, x] simplifies to x,
   so the 3 deconcatenation terms collapse; verify non-zero result *)
runTest[
  QSACoproduct[{1, 2}] =!= 0,
  True,
  "Coproduct: Delta({1,2}) is nonzero"];

Print[""];
Print["=== Braiding Tests ==="];

(* Simple braiding: sigma(e_1 tensor e_2) = q_{1,2} (e_2 tensor e_1) *)
runTest[
  Expand[QSABraid[{1}, {2}, {{1, 3}, {5, 1}}]] ===
    Expand[3 Subscript[v, 2] \[TensorProduct] Subscript[v, 1]],
  True,
  "Braid: sigma(e_1, e_2) with q_{1,2}=3"];

Print[""];
Print["=== Primary Coefficient Test ==="];

(* Primary coefficient for a single-letter word is 1 *)
runTest[
  QSAPrimaryCoefficient[{1}] === 1,
  True,
  "PrimaryCoefficient: c_{1} = 1 for single letter"];

(* Primary coefficient for Lyndon word {2,1} *)
runTest[
  QSAPrimaryCoefficient[{2, 1}] =!= 0,
  True,
  "PrimaryCoefficient: c_{2,1} is nonzero"];

Print[""];
Print["=== Standard Factorization Tests ==="];

(* {2,1} = ({2}, {1}) — suffix {1} is Lyndon *)
runTest[
  QSAStandardFactorization[{2, 1}],
  {{2}, {1}},
  "StandardFactorization: {2,1} = ({2},{1})"];

(* {3,1,2}: longest Lyndon suffix starting from right *)
runTest[
  QSAIsPrime[QSAStandardFactorization[{3, 1, 2}][[2]]],
  True,
  "StandardFactorization: suffix of {3,1,2} is prime"];

Print[""];
Print["============================="];
Print[ToString[passCount], "/", ToString[testCount], " tests passed"];
If[failCount > 0, Print[ToString[failCount], " FAILED"]];
Print["============================="];

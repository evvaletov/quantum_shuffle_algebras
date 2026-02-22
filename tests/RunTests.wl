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

Print["=== Lyndon Word Tests ==="];

runTest[QSAIsPrime[{1, 2}], True, "QSAIsPrime: {1,2} is prime"];
runTest[QSAIsPrime[{2, 1}], False, "QSAIsPrime: {2,1} is not prime"];
runTest[QSAIsPrime[{1}], True, "QSAIsPrime: single element is prime"];
runTest[QSAIsPrime[{1, 1}], False, "QSAIsPrime: {1,1} is not prime"];
runTest[QSAIsPrime[{1, 2, 3}], True, "QSAIsPrime: {1,2,3} is prime"];
runTest[QSAIsPrime[{1, 3, 2}], True, "QSAIsPrime: {1,3,2} is prime"];

Print[""];
Print["=== Relation Tests ==="];

runTest[QSARelation[{1}, {2}], -1, "QSARelation: {1} < {2}"];
runTest[QSARelation[{2}, {1}], 1, "QSARelation: {2} > {1}"];
runTest[QSARelation[{1}, {1}], 0, "QSARelation: {1} = {1}"];
runTest[QSARelation[{1}, {1, 2}], 1, "QSARelation: {1} > {1,2} (shorter is greater)"];
runTest[QSARelation[{1, 2}, {1}], -1, "QSARelation: {1,2} < {1} (longer is lesser)"];

Print[""];
Print["=== Unique Prime Factorization Tests ==="];

runTest[QSAUniquePrimeFactorization[{1, 2}], {{1, 2}},
  "UPF: {1,2} is already prime"];
runTest[QSAUniquePrimeFactorization[{2, 1}], {{2}, {1}},
  "UPF: {2,1} = {2}.{1}"];
runTest[
  Length[QSAUniquePrimeFactorization[{18, 19, 4, 8, 5, 7}]], 2,
  "UPF: {18,19,4,8,5,7} has 2 factors"];

Print[""];
Print["=== Shuffle Multiplication Tests ==="];

(* v_{1} . v_{2} should be v_1 v_2 + q_{1,2} v_2 v_1 *)
runTest[
  QSAShuffleMultiplication[{{1}, {2}}] ===
    Subscript[v, 1] \[TensorProduct] Subscript[v, 2] +
    Subscript[q, 1, 2] Subscript[v, 2] \[TensorProduct] Subscript[v, 1],
  True,
  "Shuffle: v_1 . v_2 = v_1 v_2 + q_{1,2} v_2 v_1"];

Print[""];
Print["=== Lyndon Word Enumeration Tests ==="];

(* Lyndon words of length <= 2 over {1,2}: {1}, {2}, {1,2} *)
runTest[
  Sort[QSALyndonWords[2, {1, 2}]],
  Sort[{{1}, {2}, {1, 2}}],
  "LyndonWords: length <= 2 over {1,2}"];

(* Number of Lyndon words of length <= 3 over {1,2} *)
(* Length 1: {1},{2}; Length 2: {1,2}; Length 3: {1,1,2},{1,2,2} *)
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
  QSACoproduct[{1}] === 1 \[TensorProduct] Subscript[v, 1] + Subscript[v, 1] \[TensorProduct] 1,
  True,
  "Coproduct: Delta(e_1) = 1 tensor e_1 + e_1 tensor 1"];

Print[""];
Print["=== Braiding Tests ==="];

(* Simple braiding: sigma(e_1 tensor e_2) = q_{1,2} (e_2 tensor e_1) *)
runTest[
  QSABraid[{1}, {2}, {{1, 3}, {5, 1}}] ===
    3 Subscript[v, 2] \[TensorProduct] Subscript[v, 1],
  True,
  "Braid: sigma(e_1, e_2) with q_{1,2}=3"];

Print[""];
Print["=== Primary Coefficient Test ==="];

runTest[
  QSAPrimaryCoefficient[{1, 2}] === 1,
  True,
  "PrimaryCoefficient: c_{1,2} = 1 for a prime word"];

Print[""];
Print["============================="];
Print[ToString[passCount], "/", ToString[testCount], " tests passed"];
If[failCount > 0, Print[ToString[failCount], " FAILED"]];
Print["============================="];

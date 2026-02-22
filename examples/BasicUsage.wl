(* BasicUsage.wl — Introductory examples for QuantumShuffleAlgebra *)

(* Load the package *)
Get["QuantumShuffleAlgebra`"];

(* --- 1. Lyndon word operations --- *)

Print["=== Lyndon Word Operations ==="];
Print[""];

a = {18, 19, 4, 8, 5, 7};
Print["a = ", a];

Print["Is a prime? ", QSAIsPrime[a]];

upf = QSAUniquePrimeFactorization[a];
Print["Unique prime factorization: ", upf];
Print[""];

(* --- 2. Quantum shuffle basis element --- *)

Print["=== Quantum Shuffle Basis ==="];
Print[""];

Print["X_a = ", QSAX[a]];
Print[""];

(* --- 3. Shuffle multiplication --- *)

Print["=== Shuffle Multiplication ==="];
Print[""];

b = {10, 7, 2};
c = {5, 4, 5};
Print["b = ", b];
Print["c = ", c];
Print["v_b . v_c = ", QSAShuffleMultiplication[{b, c}]];
Print[""];

(* --- 4. Express in Lyndon word basis --- *)

Print["=== Express in Lyndon Words ==="];
Print[""];

Print["Expressing v_a in terms of X_c's:"];
Print[QSAExpressInLyndonWords[a]];
Print[""];

Print["Primary coefficient alpha_a = ", QSAPrimaryCoefficient[a]];
Print[""];

(* --- 5. Lyndon word enumeration --- *)

Print["=== Lyndon Words up to Length 3 over {1,2} ==="];
Print[""];

lyndon = QSALyndonWords[3, {1, 2}];
Print[lyndon];
Print[""];

(* --- 6. Hopf algebra operations --- *)

Print["=== Hopf Algebra ==="];
Print[""];

Print["Coproduct of {1}: ", QSACoproduct[{1}]];
Print["Coproduct of {1,2}: ", QSACoproduct[{1, 2}]];
Print["Counit of {}: ", QSACounit[{}]];
Print["Counit of {1}: ", QSACounit[{1}]];
Print["Antipode of {1}: ", QSAAntipode[{1}]];
Print[""];

(* --- 7. Braiding --- *)

Print["=== Braiding ==="];
Print[""];

qMat = {{1, 2}, {3, 1}};
Print["Braiding matrix Q = ", qMat];
Print["sigma(e_1, e_2) = ", QSABraid[{1}, {2}, qMat]];

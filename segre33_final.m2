------------------------------------------------------------
------ Segre Determinant Seg(3,3)
------ Author: Lizzie Pratt
------ Date: 12.10.2025
------------------------------------------------------------


------------------------------------------------------------
----------------------- Setup ------------------------------
------------------------------------------------------------

-- Create a ring by taking the ring of the Grassmannian Gr(3,9)
-- and adding variables b_{i,j} for entries of B
-- Create a ring
-- by taking the ring of the Grassmannian Gr(3,9)
-- and adding variables b_{i,j} for entries of B

bvars = apply(9, i->b_(0,i)) | apply(9, i->b_(1,i)) | apply(9, i->b_(2,i))
pluckp = flatten flatten apply(9, i -> apply(i, j -> apply(j, k-> p_(k,j,i))))
R = QQ[bvars | pluckp];
G = Grassmannian(2,8,CoefficientRing =>QQ);
I = sub(G,R);



------------------------------------------------------------
-------------------- Helper functions -----------------------
-- Run all of these before the main computation
------------------------------------------------------------

-- Position helper
pos = (S,n) -> position(S, i -> i == n)

-- Sign computation used in Laplace expansion
sgn = (I,J,K) -> (
    S = {0,1,2,3,4,5,6,7,8};
    T = sort(S - set(I));
    I_0 + I_1 + I_2 + pos(T,J_0) + pos(T,J_1) + pos(T,J_2))


-- Laplace expansion coefficient for p_I * p_J * p_K
fullCoeff = (I,J) -> (
    K := {0,1,2,3,4,5,6,7,8} - set(I) - set(J);
    (-1)^(sgn(I,J,K))
    * b_(0,I_0)*b_(0,I_1)*b_(0,I_2)
    * b_(1,J_0)*b_(1,J_1)*b_(1,J_2)
    * b_(2,K_0)*b_(2,K_1)*b_(2,K_2)
    * p_(I_0,I_1,I_2)
    * p_(J_0,J_1,J_2)
    * p_(K_0,K_1,K_2)
)

-- test: fullCoeff({0,1,2},{3,4,5})

-- extracting coefficient of p_I*p_J*p_K
extractCoeff =(I,J,K,g) -> diff(p_(I_0,I_1,I_2),diff(p_(J_0,J_1,J_2),diff(p_(K_0,K_1,K_2),g)));

-- test: extractCoeff({0,1,2},{3,4,5},{6,7,8}, 3*p_(0,1,2)*p_(3,4,5)*p_(6,7,8))


-- Returns the vector of coefficients of bmon against all brackets
constraint = (bmon) -> apply(brackets, brack-> coefficient(bmon,brack))

-- Random constraint using a randomly selected b-monomial
randConstraint = (T) -> (
    mon = T_(random(#T));
    (constraint(leadMonomial(mon)),coefficient(leadMonomial(mon),mon)
        ))


-- Convert a polynomial in b-variables to a Plücker vector (approx via random constraints)
-- Note: this uses randomness; increase sample count if solve is unstable (currently 100)
toPlucker = (Bcoeff) -> (
    coeffTerms = terms Bcoeff;
    constraints = apply(100, i -> randConstraint(coeffTerms));
    M = matrix(apply(constraints, c -> c#0));
    N = transpose matrix{apply(constraints, c -> c#1)};
    x = solve(M,N);
    if x === null then error "solve(M,N) returned null — try increasing sample size";
    x
)
------------------------------------------------------------
--------------------- Main computation ----------------------
------------------------------------------------------------


------------------------------------------------------------
-- 1. Laplace expansion into b_{ij} and [I][J][K]
------------------------------------------------------------

S = {0,1,2,3,4,5,6,7,8};

allCoefficients =
    flatten apply(subsets(S,3), I ->
        apply(subsets(S - set(I), 3), J -> fullCoeff(I,J))
    );


f = sum(allCoefficients); 

print("Built f with #terms = " | toString (#terms f))


------------------------------------------------------------
-- 2. Straightening algorithm for Gr(3,9)
--    (takes 5–10 seconds)
------------------------------------------------------------

time h = f % I;


------------------------------------------------------------
-- 3. Define variables q_(i,j,k) indexing minors of B
------------------------------------------------------------

B = matrix{
    apply(9, i -> b_(0,i)),
    apply(9, i -> b_(1,i)),
    apply(9, i -> b_(2,i))
};

apply(subsets({0,1,2,3,4,5,6,7,8},3),I -> q_(I_0,I_1,I_2) = (det B_{I_0,I_1,I_2}));



------------------------------------------------------------
-- 4. Basis for <I><J><K> via standard Young tableaux
------------------------------------------------------------

needsPackage "SpechtModule";

ption = new Partition from {3,3,3};

tabl = standardTableaux ption;
tabl = toListOfTableaux tabl;

stdMons = apply(tabl, t -> (t^0, t^1, t^2));
print("Number of standard monomials = " | toString(#stdMons))

brackets = apply(stdMons, mon ->
    q_(mon_0#0, mon_0#1, mon_0#2)
  * q_(mon_1#0, mon_1#1, mon_1#2)
  * q_(mon_2#0, mon_2#1, mon_2#2)
);
-- DO NOT try to display brackets

------------------------------------------------------------
----------------- 5. Display utilities ----------------------
------------------------------------------------------------

-- Nice display of [I][J][K]
squareBracketString = (I,J,K) -> (
    "[" | toString(I#0) | toString(I#1) | toString(I#2) | "]"
  | "[" | toString(J#0) | toString(J#1) | toString(J#2) | "]"
  | "[" | toString(K#0) | toString(K#1) | toString(K#2) | "]"
);

-- Nice display of <I><J><K>
angleBracketString = (I,J,K) -> (
    "<" | toString(I#0) | toString(I#1) | toString(I#2) | ">"
  | "<" | toString(J#0) | toString(J#1) | toString(J#2) | ">"
  | "<" | toString(K#0) | toString(K#1) | toString(K#2) | ">"
);

-- Make the square brackets [I][J][K] into ring variables
-- So the we can multiply them by a vector of integers

bracketDisplay = apply(stdMons, t -> squareBracketString(t#0, t#1, t#2));

-- if you run this too quickly after the previous command
-- you may need to rerun 
T = QQ[bracketDisplay]; 

bracketVector = matrix{gens T};


------------------------------------------------------------
------------------ Final Coefficient List -------------------
-- Computes coefficients in <L><M><N> basis
-- Takes roughly 10–20 seconds
------------------------------------------------------------

time finalCoeffs = apply(stdMons, mon -> (
    angleBracketString(mon_0, mon_1, mon_2),
    (bracketVector * promote(toPlucker(extractCoeff(mon#0,mon#1,mon#2,h)), T))_(0,0)
));



-- Display table
scan(finalCoeffs, t ->
    print(toString(t#0) | " : " | toString(t#1))
);


------------------------------------------------------------
-- 6. Verification step
--    Compare with Segre determinant h
------------------------------------------------------------

bracketVectorInB = matrix{brackets};

-- takes about 40 seconds
time finalCoeffs =
    apply(stdMons, mon ->
        bracketVectorInB
        * toPlucker(extractCoeff(mon#0,mon#1,mon#2,h))
        * p_(mon_0#0,mon_0#1,mon_0#2)
        * p_(mon_1#0,mon_1#1,mon_1#2)
        * p_(mon_2#0,mon_2#1,mon_2#2)
    );

hConjecture = sum(finalCoeffs);

hConj = hConjecture_(0,0);

hConj == h -- should return True
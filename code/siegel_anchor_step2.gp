/* attack16 step 2: modular symbols for E=11.a3
   - algebraic symbol values phi+, phi- on paths {0, c/11}, c=1..10 and {0,oo}
   - cusp <-> torsion correspondence via Abel-Jacobi mod Lambda
   - decomposition of gamma- = a-2b as Q-linear combo of {0, c/11} symbols
*/
default(realprecision, 80);

E = ellinit([0,-1,1,0,0]);
w1 = E.omega[1]; w2 = E.omega[2];
wim = w1 - 2*w2;   /* anti-invariant period, pure imaginary 2.9176i */
Sp = msfromell(E, 1);  Mp = Sp[1]; sp = Sp[2];
Sm = msfromell(E,-1);  Mm = Sm[1]; sm = Sm[2];

ph(p) = [mseval(Mp,sp,p), mseval(Mm,sm,p)];

/* symbols {0, c/11} for c = 0..10 plus {0,oo} */
V = vector(11, c, ph([0, c/11]));
Vinf = ph([0, oo]);
print("phi {0,oo} = ", Vinf);
for(c=1,10, print("phi {0,", c, "/11} = ", V[c+1]));

/* Abel-Jacobi of cusps c/11 relative to infinity cusp:
   AJ_c = int_inf^{c/11} omega = I_c - I_inf  (mod Lambda).
   Expect: AJ_c in (1/5)Lambda/Lambda, real part (3k/5) w1, phi- part in Z. */
print("--- Abel-Jacobi of rational cusps (mod Z^2) ---");
for(c=1,10, my(j = V[c+1] - Vinf);
  print("c=", c, ":  (", j[1] - floor(j[1]), ", ", j[2] - floor(j[2]), ")   raw = ", j));

/* z(A) = 3 w1 / 5, so torsion kA <-> phi+ = 3k/5 mod 1, phi- = 0 */
print("z(A)/w1*5 = ", 5*ellpointtoz(E,[0,0])/w1);

/* candidate decomposition:
   gamma- = sum_{c=1}^5 lam_c ( {0,c/11} - {0,-c/11} )
   closed on X1(11) by construction.  Solve for period (0,1):
   sum lam_c (V[c+1] - V[12-c]) = (0, 1).  5 unknowns, 2 eqns. */
D = vector(5, c, V[c+1] - V[12-c]);
for(c=1,5, print("D_", c, " = ", D[c]));
/* relations among D_c: periods are injective on H_1, so find Q-relations */
M = mattranspose(Mat(D~));
print("ker(D) = ", matker(M));
/* particular solution: use first two components */
Msolve = Mat([D[1]~, D[2]~]~);
print("det = ", matdet(Msolve));
lam12 = Msolve^(-1) * [0,1]~;
print("lam via D1,D2 = ", lam12);
print("check period: ", lam12[1]*D[1] + lam12[2]*D[2]);

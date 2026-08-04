/* Elliptic dilogarithm numerics for (C3) proof assembly.
   E = 11.a3: y^2 + y = x^3 - x^2.  P = (0,0) generator of E(Q)=Z/5.
   Checks:
     (a) exotic relation D_E(2P) = 2 D_E(P)     (Bertin, via Brunault 3.151)
     (b) -5 D_E(P) =? 2*pi*b_11                 (our divisor computation
         (x).(y) = 5(A) - 5(2A) + winding n=1 + Bloch Thm)
     (c) D_E(P)/L(E,2) =? 11/(10 pi)            (pins Brunault (3.151))
*/
default(realprecision, 50);
E = ellinit([0,-1,1,0,0]);
w1 = E.omega[1]; w2 = E.omega[2];
w3 = w1 - w2;   /* basis with Im tau > 0 (Delta < 0: conj(w2) = w1 - w2 in lattice) */
tau = w3/w1;
printf("tau = %.30f\n", tau);
printf("check Im(tau) > 0: %d\n", imag(tau) > 0);
zP  = ellpointtoz(E, [0,0]);
z2P = ellpointtoz(E, [1,-1]);
printf("z(P)  = %.30f\n", zP);
printf("z(2P) = %.30f\n", z2P);
printf("5 z(P) mod lattice (should be ~0): %.3e\n", abs(5*zP - round(real(5*zP)) ));
write("dilog_zvals.txt", zP, " ", z2P, " ", tau);
L2 = lfun(E, 2);
b11 = 11/(4*Pi^2)*L2;
printf("L(E,2) = %.40f\n", L2);
printf("b_11   = %.40f\n", b11);
printf("target (2*pi/5)*b_11 = %.40f\n", 2*Pi/5*b11);
quit

\\ bertin_diamond.gp -- model checks + certified regulator integrals for the
\\ Bertin/Brunault symbol question (attack 13).
\\
\\ Curves:
\\   C1 cubic:  F = (X+1)(Y+1)(X+Y+1) + XY, homogenized (X+Z)(Y+Z)(X+Y+Z) + XYZ
\\   W_int:     v^2 + 4uv + v = u^3 - 2u^2 - u   (Riemann-Roch image of the C1 cubic
\\              with origin Q1=[0:1:0], via u = Y/(X+1), v = Y(X+Y)/(X+1))
\\   11.a3:     y^2 + y = x^3 - x^2  (minimal model, X1(11))
default(realprecision, 80);

W = ellinit([4,-2,1,-1,0]);
print("W_int disc = ", W.disc, "   j = ", W.j);
id = ellidentify(W);
print("ellidentify: ", id);
M = ellinit(id[1][2]);   /* [0,-1,1,0,0] = 11.a3 */
v = id[2];               /* change of variable: W = ellchangecurve(M, v) */
print("11.a3 a-invariants: ", [M.a1,M.a2,M.a3,M.a4,M.a6]);
print("change v=[u,r,s,t] = ", v);
Wc = ellchangecurve(W, v);
print("sanity: a-invariants ellchangecurve(W_int,v) = ", [Wc.a1,Wc.a2,Wc.a3,Wc.a4,Wc.a6], " disc = ", Wc.disc);

/* images of the C1-cubic points on W_int */
B  = [0,0];    /* B1 = (-1,0) on the cubic: u ~ Y/(X+1) ~ X+1 -> 0, v -> 0 */
Ad = [-1,1];   /* A1 = (0,-1):  u = -1/1, v = (-1)(0-1)/1 = 1 */
print("oncurve B: ", ellisoncurve(W,B), "  oncurve A1: ", ellisoncurve(W,Ad));
print("ord(B) on W_int: ", ellorder(W,B));
print("2B = ", ellmul(W,B,2), "  3B = ", ellmul(W,B,3), "  4B = ", ellmul(W,B,4));
print("3B == A1 image ? ", ellmul(W,B,3) == Ad);

/* transport to 11.a3 via ellchangepoint (W -> M) */
Bm  = ellchangepoint(B, v);
A1m = ellchangepoint(Ad, v);
print("on 11.a3:  B1 -> ", Bm, "   A1 -> ", A1m);
P = [0,0];
print("P multiples on 11.a3: P=", P, " 2P=", ellmul(M,P,2), " 3P=", ellmul(M,P,3), " 4P=", ellmul(M,P,4));

/* z-values: u/w1 mod 1 is isomorphism-invariant */
wM = M.omega;  wW = W.omega;
print("periods 11.a3: ", wM);
print("z(P)/w1 = ", ellpointtoz(M,P)/wM[1]);
print("z(2P)/w1 = ", ellpointtoz(M,ellmul(M,P,2))/wM[1]);
print("z(3P)/w1 = ", ellpointtoz(M,ellmul(M,P,3))/wM[1]);
print("z(4P)/w1 = ", ellpointtoz(M,ellmul(M,P,4))/wM[1]);
print("z(B1)/w1 (via W_int) = ", ellpointtoz(W,B)/wW[1]);
print("z(B1)/w1 (via 11.a3) = ", ellpointtoz(M,Bm)/wM[1]);

print("L(E,2) = ", lfun(M,2));
print("b11 = L'(E,0) = ", lfun(M,0,1));

\\ ------------------------------------------------------------------
\\ Certified integrals  int_{gamma^-} eta(f,g)  by unwrapped Riemann sums.
\\ gamma^- = anti-invariant generator: u = s*(2*w2 - w1), s in [0,1].
\\ (midpoint grid, k = 0..N-1; centered lattice representative for s > 1/2)
\\ ------------------------------------------------------------------
b11 = 0.1521471417259180494862272974786344956281;
d   = 0.1911937370843316957549544343121738161012;   /* D_E(P), dilog.py 60d */

{
int_eta_weier(N) =   /* {x,y} = Weierstrass coordinates of 11.a3 (Brunault 3.210 symbol) */
  my(L, pt, xs, ys, ax, ay, ax2, ay2, dax, day, S, k, s, z);
  L = 2*wM[2] - wM[1];
  S = 0;
  s = 0.5/N; z = s*L;
  pt = ellztopoint(M, z);
  xs = pt[1]; ys = pt[2];
  ax = arg(xs); ay = arg(ys);
  for(k = 1, N-1,
    s = (k+0.5)/N; z = s*L; if(s > 0.5, z = (s-1)*L);
    pt = ellztopoint(M, z);
    ax2 = arg(pt[1]); ay2 = arg(pt[2]);
    dax = ax2 - ax; day = ay2 - ay;
    if(dax > Pi, dax -= 2*Pi); if(dax < -Pi, dax += 2*Pi);
    if(day > Pi, day -= 2*Pi); if(day < -Pi, day += 2*Pi);
    S += log(abs(xs))*day - log(abs(ys))*dax;
    xs = pt[1]; ys = pt[2]; ax = ax2; ay = ay2;
  );
  return(S);
}

{
int_eta_c1(N) =      /* {X,Y} = coordinate functions of the C1 cubic, via W_int */
  my(L, pt, uu, vv, X, Y, xs, ys, ax, ay, ax2, ay2, dax, day, S, k, s, z);
  L = 2*wW[2] - wW[1];
  S = 0;
  s = 0.5/N; z = s*L;
  pt = ellztopoint(W, z);
  uu = pt[1]; vv = pt[2];
  X = (vv - uu^2)/(uu*(uu+1)); Y = uu*(X+1);
  xs = X; ys = Y;
  ax = arg(xs); ay = arg(ys);
  for(k = 1, N-1,
    s = (k+0.5)/N; z = s*L; if(s > 0.5, z = (s-1)*L);
    pt = ellztopoint(W, z);
    uu = pt[1]; vv = pt[2];
    X = (vv - uu^2)/(uu*(uu+1)); Y = uu*(X+1);
    ax2 = arg(X); ay2 = arg(Y);
    dax = ax2 - ax; day = ay2 - ay;
    if(dax > Pi, dax -= 2*Pi); if(dax < -Pi, dax += 2*Pi);
    if(day > Pi, day -= 2*Pi); if(day < -Pi, day += 2*Pi);
    S += log(abs(xs))*day - log(abs(ys))*dax;
    xs = X; ys = Y; ax = ax2; ay = ay2;
  );
  return(S);
}

print("\\--- certified integrals over gamma^- (anti-invariant generator) ---");
iw1 = int_eta_weier(2000);
iw2 = int_eta_weier(4000);
iw3 = int_eta_weier(8000);
print("Weierstrass {x,y} on 11.a3:");
print("  N=2000: ", iw1);
print("  N=4000: ", iw2);
print("  N=8000: ", iw3);
print("  Richardson (4*i3-i2)/3 = ", (4*iw3-iw2)/3);
print("  Brunault (3.210) value  -2 pi b11 = ", -2*Pi*b11, "  (= -5 D_E(P))");
print("  factor-1 Bloch would give -pi b11 = ", -Pi*b11,   "  (= -(5/2) D_E(P))");

ic1 = int_eta_c1(1000);
ic2 = int_eta_c1(2000);
ic3 = int_eta_c1(4000);
print("C1 cubic {X,Y}:");
print("  N=1000: ", ic1);
print("  N=2000: ", ic2);
print("  N=4000: ", ic3);
print("  target 2*pi*m(C1) = 14 pi b11 = ", 14*Pi*b11, "  (= 35 D_E(P) = 2*D_E(diamond))");
print("  factor-1 Bloch would give 7 pi b11 = ", 7*Pi*b11, "  (= (35/2) D_E(P))");
print("  i3/(14 pi b11) = ", ic3/(14*Pi*b11));

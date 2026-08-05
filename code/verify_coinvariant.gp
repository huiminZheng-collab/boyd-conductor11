\\ verify_coinvariant.gp -- Richardson verification of the coinvariant resolution.
default(realprecision, 70);
E = ellinit([0,-1,1,0,0]);
w1 = E.omega[1]; w2 = E.omega[2];
pt(u) = ellztopoint(E, u);
darg(z0, z1) = arg(z1/z0);
int_eta(W, u0, N) = {
  my(s = 0.0, P0, P1, x0, y0, x1, y1);
  P0 = pt(u0 + (0.5/N)*W); x0 = P0[1]; y0 = P0[2];
  for(j = 1, N,
    P1 = pt(u0 + ((j+0.5)/N)*W); x1 = P1[1]; y1 = P1[2];
    s += log(abs(x0))*darg(y0, y1) - log(abs(y0))*darg(x0, x1);
    P0 = P1; x0 = x1; y0 = y1;
  );
  return(s);
};
b11 = 0.15214714172591804948622729747863449562814358916422612280989;
u0 = w1/7 + w2/13;
tb = -Pi*b11;
tab = 2*Pi*b11;
print("targets: int_b = ", tb);
print("         int_{a-2b} = ", tab);
run(N) = {
  my(Ia, Ib, Iab);
  Ia = int_eta(w1, u0, N);
  Ib = int_eta(w2, u0, N);
  Iab = int_eta(w1 - 2*w2, u0, N);
  print("N = ", N);
  print("  int_a          = ", Ia);
  print("  int_b     err  = ", Ib - tb);
  print("  int_a-2b  err  = ", Iab - tab);
  return([Ia, Ib, Iab]);
};
v2 = run(2000);
v4 = run(4000);
v8 = run(8000);
v16 = run(16000);
rich2(z1, z2, z3) = (16*(4*z3 - z2)/3 - (4*z2 - z1)/3)/15;
Ra = rich2(v4[1], v8[1], v16[1]);
Rb = rich2(v4[2], v8[2], v16[2]);
Rab = rich2(v4[3], v8[3], v16[3]);
print("Richardson int_a    = ", Ra, "  (predict 0)");
print("Richardson int_b    = ", Rb);
print("  err vs -pi*b11    = ", Rb - tb);
print("Richardson int_a-2b = ", Rab);
print("  err vs 2*pi*b11   = ", Rab - tab);
print("identity: -2*int_b - int_{a-2b} = ", -2*Rb - Rab);

/* verify_coinvariant.gp -- numerical check of the coinvariant/subgroup index
   phenomenon for the eta integral of the Weierstrass symbol {x_W,y_W} on
   E: y^2+y = x^3-x^2 (11.a3).  Rev.2: centered (midpoint) discretization,
   empirical convergence order, honest residuals (no assumed error model).

   Loops are shifted by u0 = w1/7 + w2/13 to avoid zeros/poles of x_W, y_W.
   Integrand: eta = log|x| darg y - log|y| darg x, sampled as
     log|x(mid)| * [arg y(right) - arg y(left)]  etc.
   which is a centered (second-order) rule for the smooth closed loop. */
default(realprecision, 70);
E = ellinit([0,-1,1,0,0]);
w1 = E.omega[1]; w2 = E.omega[2];
pt(u) = ellztopoint(E, u);
darg(z0, z1) = arg(z1/z0);
b11 = 0.15214714172591804948622729747863449562814358916422612280989;
u0 = w1/7 + w2/13;
tb  = -Pi*b11;
tab = 2*Pi*b11;
NN = [500, 1000, 2000, 4000, 8000];

int_eta(W, N) = {
  my(s = 0.0, Pm, P0, P1);
  for(j = 0, N-1,
    Pm = pt(u0 + (j/N)*W);
    P0 = pt(u0 + ((j-0.5)/N)*W);
    P1 = pt(u0 + ((j+0.5)/N)*W);
    s += log(abs(Pm[1]))*darg(P0[2], P1[2]) - log(abs(Pm[2]))*darg(P0[1], P1[1]);
  );
  return(s);
};

err_of(j, v) = if(j==1, v[1], if(j==2, v[2]-tb, v[3]-tab));

run_all() = {
  my(V = vector(#NN));
  print("targets: int_b     = ", tb);
  print("         int_{a-2b} = ", tab);
  for(k = 1, #NN,
    V[k] = [int_eta(w1, NN[k]), int_eta(w2, NN[k]), int_eta(w1-2*w2, NN[k])];
    print("N = ", NN[k]);
    print("  int_a           = ", V[k][1]);
    print("  int_b     err   = ", err_of(2, V[k]));
    print("  int_a-2b  err   = ", err_of(3, V[k]));
  );
  return(V);
};

report_orders(V) = {
  for(j = 1, 3,
    for(k = 2, #NN,
      if(abs(err_of(j, V[k])) > 0,
        print("  comp ", j, ", N ", NN[k-1], " -> ", NN[k],
              ": p = ", log(abs(err_of(j, V[k-1])/err_of(j, V[k])))/log(2));
      );
    );
  );
};

richardson(e1, e2, p) = (2^p*e2 - e1)/(2^p - 1);

report_rich(V) = {
  my(Ra, Rb, Rab);
  Ra  = richardson(V[4][1], V[5][1], 2);
  Rb  = richardson(V[4][2], V[5][2], 2);
  Rab = richardson(V[4][3], V[5][3], 2);
  print("Richardson p=2 on finest pair N=4000/8000:");
  print("  int_a     residual vs 0        = ", Ra);
  print("  int_b     residual vs -pi*b11  = ", Rb - tb);
  print("  int_a-2b  residual vs 2*pi*b11 = ", Rab - tab);
  print("identity on extrapolants: -2*int_b - int_{a-2b} = ", -2*Rb - Rab);
};

V = run_all();
print("empirical order p = log2(err(N)/err(2N)):");
report_orders(V);
report_rich(V);

/* k=1 regulator side: z-values of torsion points on minimal model [1,-1,1,-1,0],
   tau and q for the elliptic dilogarithm (disc > 0 case). */
default(realprecision, 80);
doit() = {
  my(E, w1, w2, tau, q, zP, zT, z3P);
  E = ellinit([1,-1,1,-1,0]);
  print("disc = ", E.disc, "  conductor = ", ellglobalred(E)[1]);
  print("tors = ", elltors(E));
  w1 = E.omega[1]; w2 = E.omega[2];
  /* disc > 0: w1 real, w2 pure imaginary (negative).  Use tau = -w2/w1 > 0 imag */
  tau = -w2/w1;
  q = exp(2*Pi*I*tau);
  print("w1 = ", w1);
  print("w2 = ", w2);
  print("tau = ", tau);
  print("q = ", q);
  /* torsion: generator (0,0); 2(0,0) = ? ; (0,1) = 3(0,0) */
  print("2*(0,0) = ", ellmul(E, [0,0], 2));
  print("3*(0,0) = ", ellmul(E, [0,0], 3));
  print("4*(0,0) = ", ellmul(E, [0,0], 4));
  print("(0,1) = ?*(0,0): -(0,0) = ", ellneg(E, [0,0]));
  zP  = ellpointtoz(E, [0,0]);
  z3P = ellpointtoz(E, [0,1]);
  zT  = ellpointtoz(E, ellmul(E, [0,0], 2));
  print("z((0,0)) = ", zP, "   /w1 = ", zP/w1);
  print("z((0,1)) = ", z3P, "   /w1 = ", z3P/w1);
  print("z(2T)    = ", zT, "   /w1 = ", zT/w1, "  /w2 = ", zT/w2);
};
doit();
quit

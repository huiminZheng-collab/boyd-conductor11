default(realprecision, 60);
doit(kk) = {
  my(E, N, Em, L2, Lp0, b);
  E = ellinit(ellfromeqn(y^2 + (x^2+kk*x+1)*y + x^3));
  N = ellglobalred(E)[1];
  Em = ellminimalmodel(E);
  print("\nk=", kk, ":  conductor=", N, "  j=", E.j);
  print("   minimal model ainvs = ", [Em.a1, Em.a2, Em.a3, Em.a4, Em.a6]);
  L2 = lfun(Em, 2);
  Lp0 = lfun(Em, 0, 1);
  b = N/(4*Pi^2)*L2;
  print("   L(E,2)  = ", L2);
  print("   L'(E,0) = ", Lp0);
  print("   N/(4Pi^2)*L(E,2) = ", b, "   |diff| = ", abs(Lp0-b));
};
doit(-2); doit(-1); doit(0); doit(1); doit(2);
quit

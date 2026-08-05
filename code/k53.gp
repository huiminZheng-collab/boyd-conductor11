default(realprecision, 60);
{
  my(E, Em, N, om, w1, w2, L2, b);
  E = ellinit(ellfromeqn(y^2 + (x^2-x+1)*y + x^3));
  N = ellglobalred(E)[1];
  Em = ellminimalmodel(E);
  print("conductor = ", N, "   minimal ainvs = ", [Em.a1,Em.a2,Em.a3,Em.a4,Em.a6]);
  print("disc = ", Em.disc, "  (Delta<0 => H_1^- rank 1)");
  om = Em.omega;
  w1 = om[1]; w2 = om[2];
  print("w1 = ", w1);
  print("w2 = ", w2);
  print("w_anti = 2i*Im(w2) = ", 2*I*imag(w2));
  L2 = lfun(Em, 2);
  b = lfun(Em, 0, 1);
  print("L(E,2) = ", L2);
  print("b_53 = L'(E,0) = ", b);
  print("N/(4Pi^2)*L2  = ", N/(4*Pi^2)*L2, "   diff = ", abs(b - N/(4*Pi^2)*L2));
}
quit

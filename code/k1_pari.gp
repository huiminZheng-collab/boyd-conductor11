/* k=1 (conductor 17) basic identification, mirrors n1_certify.gp + dilog.gp */
default(realprecision, 80);
doit() = {
  my(E, Emin, gr, om, w1, w2, wanti, A, Q, zA, z2A, w3, tau, L2, b17);
  E = ellinit(ellfromeqn(y^2 + (x^2+x+1)*y + x^3));
  print("model from ellfromeqn: [a1,a2,a3,a4,a6] = ", E[1..5]);
  print("disc = ", E.disc);
  print("j    = ", E.j);
  Emin = ellminimalmodel(E);
  print("minimal model: ", Emin[1..5]);
  gr = ellglobalred(Emin);
  print("conductor = ", gr[1]);
  print("tors = ", elltors(Emin));
  /* map (0,0) of original model: find via ellfromeqn structure?
     Instead: identify 17.a curve, then use torsion points on minimal model. */
  om = Emin.omega; w1 = om[1]; w2 = om[2];
  wanti = 2*I*imag(w2);
  print("w1 = ", w1);
  print("w2 = ", w2);
  print("w_anti = ", wanti);
  print("w_anti (40 digits) = ", precision(wanti,40));
  /* torsion points on minimal model */
  my(T);
  T = elltors(Emin);
  print("torsion structure: ", T);
  /* L-values */
  L2 = lfun(Emin, 2);
  b17 = 17/(4*Pi^2)*L2;
  print("L(E,2) = ", L2);
  print("b_17 = 17/(4 pi^2) L(E,2) = ", b17);
  print("root number check, L(E,1) = ", lfun(Emin, 1));
  /* q-parameter with Delta<0 convention w3 = w1 - w2 */
  w3 = w1 - w2;
  tau = w3/w1;
  print("tau = ", tau, "  Im>0: ", imag(tau) > 0);
};
doit();
quit

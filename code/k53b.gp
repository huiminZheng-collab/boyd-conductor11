default(realprecision, 60);
{
  my(E0, Em, om0);
  E0 = ellinit(ellfromeqn(y^2 + (x^2-x+1)*y + x^3));   /* cubic model k=-1 */
  Em = ellminimalmodel(E0);
  om0 = E0.omega;
  print("cubic model periods: w1' = ", om0[1]);
  print("                     w2' = ", om0[2]);
  print("w_anti' = 2i*Im(w2') = ", 2*I*imag(om0[2]));
  print("model constant kappa = w1'/w1_min = ", om0[1]/Em.omega[1]);
  print("P(L1)/w_anti' = ", -1.69246409695100638*I/(2*I*imag(om0[2])));
  print("torsion = ", elltors(Em));
  print("rank = ", ellrank(Em));
  print("(0,0) on curve? ", ellisoncurve(E0,[0,0]));
  print("order of (0,0): ", ellorder(E0,[0,0]));
}
quit

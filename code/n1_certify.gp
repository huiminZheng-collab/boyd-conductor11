default(realprecision, 80);
{
  my(E, P, Pb, z, zb, zcj, mu, om, w1, w2, wanti);
  E = ellinit(ellfromeqn(y^2 + (x^2+1)*y + x^3));   /* S_0 cubic model */
  om = E.omega; w1 = om[1]; w2 = om[2];
  wanti = 2*I*imag(w2);
  print("w1 = ", w1);
  print("w2 = ", w2);
  print("w_anti = ", wanti);
  P  = [I, exp(I*Pi/4)];                 /* crossing point (i, e^{i pi/4}) */
  Pb = [ -I, exp(-I*Pi/4)];              /* conjugate */
  print("P on curve:  ", ellisoncurve(E, P));
  print("Pb on curve: ", ellisoncurve(E, Pb));
  print("P neg = -P?  ", ellneg(E,P) == [I, -exp(I*Pi/4)]);
  z  = ellpointtoz(E, P);
  zb = ellpointtoz(E, Pb);
  zcj = conj(z);
  print("z(P)  = ", z);
  print("z(Pb) = ", zb);
  print("conj(z(P)) = ", zcj);
  /* mu = zb - conj(z) reduced mod lattice: expect (small integer combo) */
  mu = zb - zcj;
  print("mu = zb - conj(z) = ", mu);
  print("mu/w1 = ", mu/w1);
  print("(mu - round(real(mu/w1))*w1)/w2 = ", (mu - round(real(mu/w1))*w1)/w2);
  print("Im z(P) = ", imag(z));
  print("8i Im z(P) = ", 8*I*imag(z));
  print("ratio 8 Im z(P) / (2 Im w2) = ", 8*imag(z)/(2*imag(w2)));
  print("ratio 2 Im z(P) / Im w2     = ", 2*imag(z)/imag(w2));
  print("Re z(P)/w1 = ", real(z)/w1);
}
quit

/* Check ellfromeqn output for k=0 vs k=1, and locate torsion/corner points */
default(realprecision, 60);
doit() = {
  my(E0, E1, P, Pt);
  E0 = ellinit(ellfromeqn(y^2 + (x^2+1)*y + x^3));
  print("k=0 model: ", E0[1..5], " disc = ", E0.disc);
  print("k=0: [I, exp(I pi/4)] on curve: ", ellisoncurve(E0, [I, exp(I*Pi/4)]));
  print("k=0: [0,0] on curve: ", ellisoncurve(E0, [0,0]));
  if(ellisoncurve(E0,[0,0]), print("k=0: ord(0,0) = ", ellorder(E0, [0,0])));
  E1 = ellinit(ellfromeqn(y^2 + (x^2+x+1)*y + x^3));
  print("k=1 model: ", E1[1..5], " disc = ", E1.disc);
  print("k=1: [0,0] on curve: ", ellisoncurve(E1, [0,0]));
  if(ellisoncurve(E1,[0,0]), print("k=1: ord(0,0) = ", ellorder(E1, [0,0])));
  /* corner point on S_1: x = exp(2 pi i/3) = w, y = i  (y^2 = -1) */
  my(w);
  w = exp(2*I*Pi/3);
  print("k=1: (w, i) on S_1 model? ",
        abs((I)^2 + (w^2+w+1)*I + w^3) < 1e-40);
  print("k=1: [w, I] on ellfromeqn curve: ", ellisoncurve(E1, [w, I]));
  /* torsion points on E1 model */
  my(v);
  v = elltors(E1);
  print("k=1 tors on E1: ", v);
  if(#v[3] > 0, P = v[3][1]; print("generator P = ", P);
     print("2P = ", ellmul(E1, P, 2));
     print("3P = ", ellmul(E1, P, 3));
     print("4P = ", ellmul(E1, P, 4)));
};
doit();
quit

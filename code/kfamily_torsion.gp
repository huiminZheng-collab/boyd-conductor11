default(realprecision, 60);
doit(kk) = {
  my(E0, Em, N, tors, r, o);
  E0 = ellinit(ellfromeqn(y^2 + (x^2+kk*x+1)*y + x^3));
  N = ellglobalred(E0)[1];
  Em = ellminimalmodel(E0);
  tors = elltors(Em)[1];
  r = ellrank(Em)[1];
  o = ellorder(E0, [0,0]);
  print("k=", kk, "  N=", N, "  tors=Z/", tors, "  rank=", r,
        "  ord(0,0)=", if(o==0, "INFINITE", o));
};
doit(-3); doit(-2); doit(-1); doit(0); doit(1); doit(2); doit(3);
quit

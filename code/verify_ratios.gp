default(realprecision, 80);
\\ ntilde(k)/m(k) values from ntilde_family.py (mpmath, 50 dps) -- hardcoded 45 digits
ntm3 = 0.970767075929641590831458;
ntm2 = 0.721282593073082217909532;
ntm1 = -0.46493780021441118890735356443424857552297430198564;
nt0  = -0.15214714172591804948622729747863449562814358916423;
nt1  = 0.29935558688291539005379974769003145062973595594494;
nt2  = 0.71524093225499529999788513842759388772993771089174;
nt3  = 1.07149673434035160705377;
mm2  = 0.721282593073082217909532;
mm1  = 0.50926262288824211511625928008401754294373479778084;
m0   = 0.405602955915010403908189961285;
m1   = 0.506562114988019791399529;
m2   = 0.71524093225499529999788513842759388772993771089174;
m3   = 1.07149673434035160705377;

rat(v, b) = { my(r); if(b==0, return("  b=0")); r = v/b;
  return(concat(Strprintf("%.30f", r), concat("   bestappr: ", Str(bestappr(r, 10^6))))); };

doit(kk, nt, mm) = {
  my(E, N, Em, b);
  E = ellinit(ellfromeqn(y^2 + (x^2+kk*x+1)*y + x^3));
  N = ellglobalred(E)[1];
  Em = ellminimalmodel(E);
  b = abs(lfun(Em, 0, 1));      \\ |L'(E,0)|, sign tells root number
  print("\nk=", kk, "  N=", N, "  j=", E.j, "  L'(E,0)=", lfun(Em,0,1));
  print("   |nt|/b: ", rat(abs(nt), b));
  print("   m/b   : ", rat(mm, b));
  print("   lindep[nt,b]: ", lindep([nt, b]));
  print("   lindep[m ,b]: ", lindep([mm, b]));
};
doit(-3, ntm3, m3*0+0.970767075929641590831458);
doit(-2, ntm2, mm2);
doit(-1, ntm1, mm1);
doit(0, nt0, m0);
doit(1, nt1, m1);
doit(2, nt2, m2);
doit(3, nt3, m3);
quit

/* attack16 step 1: basics on E: y^2+y=x^3-x^2 (11.a3)
   - torsion points, periods, anti-invariant cycle
   - independent numerical verification of div(x), div(y) for the
     Boyd functions transported to E via kappa_exact.py map:
       x_E = (X+Y)/(X-1)
       y_E = (X-1) - (X+Y)(Y+1)/(X-1)^2
*/
default(realprecision, 80);

E = ellinit([0,-1,1,0,0]);
print("disc = ", E.disc, "  j = ", E.j);
print("omega = ", E.omega);
w1 = E.omega[1]; w2 = E.omega[2];
print("w1 - 2 w2 = ", w1 - 2*w2);
T = elltors(E);
print("tors = ", T);
A = [0,0];
for(k=2,4, Ak = ellmul(E,A,k); print(k,"A = ", Ak));
print("5A = ", ellmul(E,A,5));

/* order of rational function f(X,Y) at point P (affine, 2Y+1 != 0) */
ordat(E, f, P) = {
  my(t = 1e-8, xs, ys, y0 = P[2], v1, v2);
  xs = P[1] + t;
  ys = ellordinate(E, xs);
  /* pick branch near y0 */
  if(abs(ys[1]-y0) < abs(ys[2]-y0), v1 = subst(subst(f,X,xs),Y,ys[1]),
                                  v1 = subst(subst(f,X,xs),Y,ys[2]));
  xs = P[1] + 10*t;
  ys = ellordinate(E, xs);
  if(abs(ys[1]-y0) < abs(ys[2]-y0), v2 = subst(subst(f,X,xs),Y,ys[1]),
                                  v2 = subst(subst(f,X,xs),Y,ys[2]));
  return(round(log(abs(v1/v2))/log(1/10)));
}

/* order at O: use v(x) = -2 */
ordatO(E, f) = {
  my(xs = 1e20, ys, v1, v2);
  ys = ellordinate(E, xs); v1 = subst(subst(f,X,xs),Y,ys[1]);
  xs = 10*xs;
  ys = ellordinate(E, xs); v2 = subst(subst(f,X,xs),Y,ys[1]);
  return(round(2*log(abs(v1/v2))/log(10)));
}

fxE = (X+Y)/(X-1);
fyE = ((X-1)^3 - (X+Y)*(Y+1))/(X-1)^2;

print("--- div(x_E), expect [A]+[2A]-[O]-[3A] ---");
print("ord at A  = ", ordat(E,fxE,[0,0]));
print("ord at 2A = ", ordat(E,fxE,ellmul(E,A,2)));
print("ord at 3A = ", ordat(E,fxE,ellmul(E,A,3)));
print("ord at 4A = ", ordat(E,fxE,ellmul(E,A,4)));
print("ord at O  = ", ordatO(E,fxE));

print("--- div(y_E), expect 3[2A]-2[3A]-[O] ---");
print("ord at A  = ", ordat(E,fyE,[0,0]));
print("ord at 2A = ", ordat(E,fyE,ellmul(E,A,2)));
print("ord at 3A = ", ordat(E,fyE,ellmul(E,A,3)));
print("ord at 4A = ", ordat(E,fyE,ellmul(E,A,4)));
print("ord at O  = ", ordatO(E,fyE));

/* z-coordinate of A */
print("z(A) = ", ellpointtoz(E, A));
print("z(A)/(w1/5) = ", ellpointtoz(E,A)/(w1/5));

/* b11 reference */
print("L'(E,0) = ", lfun(E,0,1));
print("2 Pi L'(E,0) = ", 2*Pi*lfun(E,0,1));
print("L(E,1) = ", lfun(E,1));

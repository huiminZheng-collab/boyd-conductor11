default(realprecision, 60);
E = ellinit([0,-1,1,0,0]);  /* 11a3 = X_1(11) */
w = E.omega;
print("w1 = ", w[1]); print("w2 = ", w[2]);
zA = ellpointtoz(E, [0,0]);
print("z(A) = ", zA, "  3w1/5 = ", 3*w[1]/5, "  diff mod lattice check: ", (zA - 3*w[1]/5));
NC = 500;
an = ellan(E, NC);
z11 = exp(2*Pi*I/11);
B2(x) = x^2 - x + 1/6;
gunit(a, b, tau) = {
  my(at = a % 11);
  exp(Pi*I*tau*B2(at/11))
  * prod(n=0, 120, 1 - z11^b*exp(2*Pi*I*tau*(n+at/11)))
  * prod(n=1, 120, 1 - z11^(-b)*exp(2*Pi*I*tau*(n-at/11)));
};
logG(a, tau) = sum(b=0, 10, log(gunit(a, b, tau)));
zval(tau) = sum(n=1, NC, an[n]*exp(2*Pi*I*n*tau)/n);
{
  foreach([I, 2*I, I+1/2, (3*I+1)/2], t,
    my(P = ellztopoint(E, zval(t)), X = P[1], Y = P[2], xB, yB, rx, ry);
    xB = (X+Y)/(X-1);
    yB = ((X-1)^3 - (X+Y)*(Y+1))/(X-1)^2;
    rx = log(xB) - (logG(4,t)+logG(5,t)-2*logG(2,t));
    ry = log(yB) - (logG(1,t)-3*logG(2,t)-logG(3,t)+3*logG(5,t));
    print("tau = ", t);
    print("  logCx = ", rx);
    print("  logCy = ", ry));
}

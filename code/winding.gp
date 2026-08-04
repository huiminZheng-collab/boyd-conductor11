default(realprecision, 60);
E = ellinit([0,-1,1,0,0]);
w1 = E.omega[1]; w2 = E.omega[2];
print("w1 = ", w1);
print("w2 = ", w2);
print("tau = ", w2/w1);
Om_re = 6.346046521397767108443973083772720814022;
Om_im = -2.917633233876990458661779225807326303805*I;
I_loop = -0.4744703643674631282653731076758462739531*I;
kappa = w1 / Om_re;
Imin = kappa * I_loop;
print("kappa = ", kappa);
print("I_loop on minimal model = ", Imin);
w_anti = 2*I*imag(w2);
print("w_anti = ", w_anti);
print("ratio = ", Imin / w_anti);
print("bestappr(ratio, 10^6) = ", bestappr(Imin / w_anti, 10^6));
quit

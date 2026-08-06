# DIRECT numerical integration of int eta(u,v) along Manin symbols -- no Brunault Lambda.
# Decides the pi-normalization question:
#   (a) cond-14 Brunault example: paper says int_{-2/7}^{2/7} eta(u,v) = Lambda(4f14,0) = 4 L'(E14,0)
#       = 0.9099248920494044...  (pi-times alternative: 4*pi*L' = 2.85772...)
#   (b) cond-11 main: int_{gamma^-} eta(x,y);  Siegel/Lambda chain says -2 b11 = -0.30429...
#       repo regulator chain says +/-2 pi b11 = +/-0.95597...
# eta(f,g) = log|f| darg g - log|g| darg f,  darg z = d(arg z).
# Path: each Manin symbol (alpha) = {alpha*0, alpha*iinf} -> tau = i*t, t in (0, inf).
# t>=1: direct series; t in (0,1]: Lemma 4 S-transform g_{a,b}(i/t) = w g_{b,-a}(i/(1/t)).
from mpmath import mp, mpf, mpc, pi, exp, log, sin, cos, nstr, quad

mp.dps = 50
I = mpc(0, 1)

def B2(x): return x*x - x + mpf(1)/6

class Units:
    # function F(tau) = const * prod_j g_{a_j,b_j}(tau)^{n_j}; only log|F|, darg F matter
    def __init__(self, N, terms):   # terms = [(a,b,n), ...]
        self.N = N
        self.terms = [(a % N, b % N, n) for (a, b, n) in terms]
        self.z = exp(2*pi*I/N)
    def transformed(self, al):      # g o alpha: (a,b) -> (a,b) alpha  (mod N)
        (r, s, t, u) = al
        return Units(self.N, [((a*r + b*t) % self.N, (a*s + b*u) % self.N, n)
                              for (a, b, n) in self.terms])
    def S(self):                    # Lemma 4: g_{a,b}(i/t) = w g_{b,-a}(i/(1/t))
        return Units(self.N, [(b, (-a) % self.N, n) for (a, b, n) in self.terms])
    def logabs_darg(self, u):
        # returns (log|F(i*u)|, d/d u arg F(i*u)) for u >= 1
        N = self.N; z = self.z
        lg = mpf(0); dg = mpf(0)
        q = exp(-2*pi*u)
        for (a, b, n_) in self.terms:
            at = a % N
            lg += n_ * (-pi * B2(mpf(at)/N) * u)
            # product 1: n >= 0: (1 - z^b q^{n + at/N}); product 2: n >= 1: (1 - z^{-b} q^{n - at/N})
            for (start, bb, off) in ((0, b, mpf(at)/N), (1, (-b) % N, -mpf(at)/N)):
                n = start
                while True:
                    x = n + off
                    if x <= 0:
                        n += 1; continue
                    qx = q**x if x != 0 else mpf(1)
                    if x > 0 and qx < mpf(10)**(-55):
                        break
                    w = z**bb * qx
                    lg += n_ * log(abs(1 - w))
                    c = 2*pi*x
                    dg += n_ * ((w * c) / (1 - w)).imag
                    n += 1
        return lg, dg

def eta_val(F, G, u):
    lf, df = F.logabs_darg(u)
    lg, dg = G.logabs_darg(u)
    return lf*dg - lg*df

def symbol_int(F, G, al):
    Fa, Ga = F.transformed(al), G.transformed(al)
    Fs, Gs = Fa.S(), Ga.S()
    # int_1^inf eta(Fa,Ga)(iu) du + int_0^1 eta(Fa,Ga)(it) dt
    # second = int_1^inf eta(Fs,Gs)(iu) du   (S-transform, verified substitution)
    v1 = quad(lambda u: eta_val(Fa, Ga, u), [1, 8, 40, 160])
    v2 = quad(lambda u: eta_val(Fs, Gs, u), [1, 8, 40, 160])
    return v1 + v2

print("="*70)
print("(a) cond-14: u = g5 g6/(g1 g2), v = -g3 g5 g6^2/(g1^2 g2 g4)")
u14 = Units(14, [(0,5,1),(0,6,1),(0,1,-1),(0,2,-1)])
v14 = Units(14, [(0,3,1),(0,5,1),(0,6,2),(0,1,-2),(0,2,-1),(0,4,-1)])
sym14 = [(-1, (2,1,7,4)), (-1, (1,0,4,1)), (1, (1,0,-4,1)), (1, (-2,1,7,-4))]
tot14 = mpf(0)
for sgn, al in sym14:
    v = symbol_int(u14, v14, al)
    print("  symbol (%+d) %s: %s" % (sgn, al, nstr(v, 40)), flush=True)
    tot14 += sgn*v
print("int_{-2/7}^{2/7} eta(u,v) =", nstr(tot14, 45))
print("  paper: Lambda(4f14,0) = 4 L'(E14,0) = 0.90992489204940443157992925809053")
print("  pi*paper                                = 2.85772363626847739...")

print("="*70)
print("(b) cond-11: x = -G4 G5/G2^2, y = G1 G5^3/(G2^3 G3), gamma^- 7 symbols")
x11 = Units(11, [(2,b,-2) for b in range(11)] + [(4,b,1) for b in range(11)] + [(5,b,1) for b in range(11)])
y11 = Units(11, [(1,b,1) for b in range(11)] + [(2,b,-3) for b in range(11)]
               + [(3,b,-1) for b in range(11)] + [(5,b,3) for b in range(11)])
sym11 = [(1, (1,0,3,1)), (-1, (1,1,3,4)), (1, (3,1,11,4)),
         (-1, (1,0,1,1)), (1, (1,2,1,3)), (-1, (3,2,4,3)), (1, (3,8,4,11))]
tot11 = mpf(0)
for sgn, al in sym11:
    v = symbol_int(x11, y11, al)
    print("  symbol (%+d) %s: %s" % (sgn, al, nstr(v, 40)), flush=True)
    tot11 += sgn*v
print("int_{gamma^-} eta(x,y) =", nstr(tot11, 45))
print("  -2 b11      = -0.30429428345183609897245459495727")
print("  +2 pi b11   = +0.95596868542165847877477217156")
print("  -2 pi b11   = -0.95596868542165847877477217156")

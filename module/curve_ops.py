import curve25519_python as cp
from utils import E, bump, rand

B = b'\x58' + b'\x66'*31
ONE = b'\x01' + b'\0'*31

def P(A, B):
    return E(cp.point_addition(A, B))

def Pm(s, P_):
    bump()
    return E(cp.point_multiply(s, P_))

def Sm(s):
    bump()
    return E(cp.scalar_multiply(s))

def S2(a, b):
    return E(cp.scalar_multiply_scalar(a, b))

def Inv(z):
    return E(cp.scalar_inverse(z))

def Sub(a, b):
    return E(cp.scalar_subtraction(a, b))

H1 = Pm(b'\x02'+b'\0'*31, B)
H2 = Pm(b'\x03'+b'\0'*31, B)

def cm(r, s):
    return P(Pm(r, H1), Pm(s, H2))

def rho(h, r, s):
    return P(Pm(r, h), Pm(s, B))
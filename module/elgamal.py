from curve_ops import Pm, P, Inv, B, rand

def elgamal_encrypt(m, pk, sk=None):
    """Encrypt m (point) under pk (public key point) on Curve25519. If sk is provided, include decryption key."""
    k = rand()  # Random scalar
    c1 = Pm(k, B)  # k * G
    c2 = P(m, Pm(k, pk))  # m + k * pk
    if sk:
        return (c1, c2, sk)  # Include sk for dealer
    return (c1, c2)

def elgamal_decrypt(c1, c2, sk):
    """Decrypt (c1, c2) using sk to recover m."""
    temp = Pm(sk, c1)  # sk * c1 = sk * k * G
    return P(c2, Pm(Inv(sk), temp))  # c2 - sk * c1 = m
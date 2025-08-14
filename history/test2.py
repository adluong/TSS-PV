import secrets
import hashlib
import curve25519_python as cp

IDENTITY_POINT = b'\x01' + b'\x00' * 31

def ensure_bytes(data):
    if isinstance(data, list):
        return bytes(data)
    elif isinstance(data, bytes):
        return data
    else:
        raise TypeError(f"Expected bytes or list, got {type(data)}")

def random_scalar():
    return secrets.token_bytes(32)

def random_point():
    scalar = random_scalar()
    return ensure_bytes(cp.scalar_multiply(scalar))

# Schnorr ZKP for f = h^r + w^s
def zkpok_proof(h: bytes, w: bytes, f: bytes, r: bytes, s: bytes):
    y = random_scalar()
    z = random_scalar()
    h_y = ensure_bytes(cp.point_multiply(y, h))
    w_z = ensure_bytes(cp.point_multiply(z, w))
    commitment = ensure_bytes(cp.point_addition(h_y, w_z))
    hash_data = f + commitment
    hash_value = hashlib.sha256(hash_data).digest()[:32]
    hash_r = ensure_bytes(cp.scalar_multiply_scalar(hash_value, r))
    hash_s = ensure_bytes(cp.scalar_multiply_scalar(hash_value, s))
    u = ensure_bytes(cp.scalar_subtraction(y, hash_r))
    v = ensure_bytes(cp.scalar_subtraction(z, hash_s))
    return (hash_value, u, v, f)

def zkpok_verify(h: bytes, w: bytes, proof):
    hash_value, u, v, f = proof
    f_hash = ensure_bytes(cp.point_multiply(hash_value, f))
    h_u = ensure_bytes(cp.point_multiply(u, h))
    w_v = ensure_bytes(cp.point_multiply(v, w))
    med = ensure_bytes(cp.point_addition(f_hash, h_u))
    M = ensure_bytes(cp.point_addition(med, w_v))
    hash_data = f + M
    hash_value_prime = hashlib.sha256(hash_data).digest()[:32]
    return hash_value_prime == hash_value

def generate_shares_with_zkp(secret_scalar: bytes, n: int, k: int):
    coefficients = [secret_scalar] + [secrets.token_bytes(32) for _ in range(k - 1)]
    x_values = [int.to_bytes(i + 1, 32, 'little') for i in range(n)]
    shares = []
    proofs = []
    for x in x_values:
        # Classic Shamir polynomial evaluation in the field
        y = coefficients[0]
        x_pow = x
        for coef in coefficients[1:]:
            term = ensure_bytes(cp.scalar_multiply_scalar(coef, x_pow))
            y = ensure_bytes(cp.scalar_addition(y, term))
            x_pow = ensure_bytes(cp.scalar_multiply_scalar(x_pow, x))
        # Shamir share as a point
        h = ensure_bytes(cp.scalar_multiply(y))
        # Blinding generator
        w = random_point()
        # Randomizer scalars
        r = random_scalar()
        s = random_scalar()
        # f = h^r + w^s
        hr = ensure_bytes(cp.point_multiply(r, h))
        ws = ensure_bytes(cp.point_multiply(s, w))
        f = ensure_bytes(cp.point_addition(hr, ws))
        # ZKP of knowledge of r, s s.t. f = h^r + w^s
        proof = zkpok_proof(h, w, f, r, s)
        shares.append((x, f, h, w, r, s))
        proofs.append(proof)
    return shares, proofs

def reconstruct_shares(shares, proofs, k):
    # Use the first k valid ZKP shares
    valid_shares = []
    for share, proof in zip(shares, proofs):
        x, f, h, w, r, s = share
        if zkpok_verify(h, w, proof):
            valid_shares.append((x, f))
        if len(valid_shares) == k:
            break
    if len(valid_shares) < k:
        raise ValueError("Not enough valid shares to reconstruct!")
    return reconstruct_original(valid_shares)

def reconstruct_original(shares):
    k = len(shares)
    total_yj_delta = IDENTITY_POINT
    for j in range(k):
        xj_bytes, yj_bytes = shares[j]
        delta_bytes = b'\x01' + b'\x00' * 31
        for m in range(k):
            if m != j:
                xm_bytes, _ = shares[m]
                xm_minus_xj = ensure_bytes(cp.scalar_subtraction(xm_bytes, xj_bytes))
                inv_xm_minus_xj = ensure_bytes(cp.scalar_inverse(xm_minus_xj))
                term = ensure_bytes(cp.scalar_multiply_scalar(xm_bytes, inv_xm_minus_xj))
                delta_bytes = ensure_bytes(cp.scalar_multiply_scalar(delta_bytes, term))
        scaled_yj = ensure_bytes(cp.point_multiply(delta_bytes, yj_bytes))
        total_yj_delta = ensure_bytes(cp.point_addition(total_yj_delta, scaled_yj))
    return total_yj_delta

def test_sharing_with_zkp():
    n, k = 32, 22
    secret_scalar = random_scalar()
    original_secret = ensure_bytes(cp.scalar_multiply(secret_scalar))
    print(f"Original:      {original_secret.hex()}")
    print(f"Generating {n} shares (with ZKP as in comp.py)...")
    shares, proofs = generate_shares_with_zkp(secret_scalar, n, k)
    print("Reconstructing from first k valid shares...")
    reconstructed_secret = reconstruct_shares(shares, proofs, k)
    print(f"Reconstructed: {reconstructed_secret.hex()}")
    if reconstructed_secret == original_secret:
        print("Test passed: Reconstructed secret matches the original.")
    else:
        print("Test failed: Reconstructed secret does not match the original.")

if __name__ == "__main__":
    test_sharing_with_zkp()

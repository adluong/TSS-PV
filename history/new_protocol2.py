import secrets
import hashlib
import curve25519_python as cp
import time
import psutil
import os
IDENTITY_POINT = b'\x01' + b'\x00' * 31

GENERATOR_POINT = (
    b'\x58\x66\x66\x66\x66\x66\x66\x66'
    b'\x66\x66\x66\x66\x66\x66\x66\x66'
    b'\x66\x66\x66\x66\x66\x66\x66\x66'
    b'\x66\x66\x66\x66\x66\x66\x66\x66'
)
n = 32
k = 2*n//3 + 1
# k = n//3 + 1
# Helper functions
def ensure_bytes(data):
    """
    Convert data to bytes if it’s a list, or return it unchanged if it’s already bytes.
    Raise an error if it’s neither.
    """
    if isinstance(data, list):
        return bytes(data)
    elif isinstance(data, bytes):
        return data
    else:
        raise TypeError(f"Expected bytes or list, got {type(data)}")

# public fixed bases h1 , h2  (here: 2·B and 3·B so they are independent of g=B)
H1_POINT = ensure_bytes(cp.point_multiply(b'\x02' + b'\x00' * 31, GENERATOR_POINT))
H2_POINT = ensure_bytes(cp.point_multiply(b'\x03' + b'\x00' * 31, GENERATOR_POINT))


def random_scalar():
    return secrets.token_bytes(32)

def random_point():
    scalar = random_scalar()
    return ensure_bytes(cp.scalar_multiply(scalar))

# Core functions with ensured bytes
def gen_glist(t: int) -> list[bytes]:
    return [ensure_bytes(cp.scalar_multiply(random_scalar())) for _ in range(t)]

def compute_hp_no_mask(x: bytes, glist: list[bytes]) -> bytes:
    hp = IDENTITY_POINT
    x_pow = x
    for g in glist:
        p_k = ensure_bytes(cp.point_multiply(x_pow, g))
        hp = ensure_bytes(cp.point_addition(hp, p_k))
        x_pow = ensure_bytes(cp.scalar_multiply_scalar(x_pow, x))
    return hp

def compute_fp(hlist, wlist, r: bytes, s: bytes) -> list[bytes]:
    if len(hlist) != len(wlist):
        raise ValueError("hlist and wlist must have the same length")
    flist = []
    for i in range(len(hlist)):
        hr = ensure_bytes(cp.point_multiply(r, hlist[i]))
        ws = ensure_bytes(cp.point_multiply(s, wlist[i]))
        fp_i = ensure_bytes(cp.point_addition(hr, ws))
        flist.append(fp_i)
        # print(f"flist: {flist[i].hex()}")
    return flist

def reconstruct(shares: list[tuple[bytes, bytes]]) -> bytes:
    k = len(shares)
    total_yj_delta = IDENTITY_POINT
    for j in range(k):
        t = 0
        xj_bytes, yj_bytes = shares[j]
        delta_bytes = IDENTITY_POINT
        start = time.perf_counter()
        for m in range(k):
            if m != j:
                xm_bytes, _ = shares[m]
                xm_minus_xj = ensure_bytes(cp.scalar_subtraction(xm_bytes, xj_bytes))
                inv_xm_minus_xj = ensure_bytes(cp.scalar_inverse(xm_minus_xj))
                term = ensure_bytes(cp.scalar_multiply_scalar(xm_bytes, inv_xm_minus_xj))
                delta_bytes = ensure_bytes(cp.scalar_multiply_scalar(delta_bytes, term))
        stop = time.perf_counter()
        t = t + (stop - start)
        # print(f"========> time to compute xj/xi-xj etc. {t:.3f}")
        scaled_yj = ensure_bytes(cp.point_multiply(delta_bytes, yj_bytes))
        total_yj_delta = ensure_bytes(cp.point_addition(total_yj_delta, scaled_yj))
    print(f" - There are {k} muls in reconstruct\n")   
    return total_yj_delta
        ##########################################################

def zkpok_proof(h_values: list[bytes], w_values: list[bytes], f_values: list[bytes], r: bytes, s: bytes) -> tuple[bytes, bytes, bytes, list[bytes]]:
    # Sample random scalars y and z as 32-byte objects
    y = secrets.token_bytes(32)
    z = secrets.token_bytes(32)

    # Compute commitments for each index
    commitments = []
    for i in range(len(h_values)):
        # h_i^y = point_multiply(y, h_values[i])
        h_i_y = ensure_bytes(cp.point_multiply(y, h_values[i]))
        # w_i^z = point_multiply(z, w_values[i])
        w_i_z = ensure_bytes(cp.point_multiply(z, w_values[i]))
        # commitment = h_i^y + w_i^z
        commitment = ensure_bytes(cp.point_addition(h_i_y, w_i_z))
        commitments.append(commitment)

    # Create hash input by concatenating f_values and commitments as bytes
    hash_data = b"".join(f_values + commitments)
    # Compute SHA-256 hash and truncate to 32 bytes for scalar compatibility
    hash_value = hashlib.sha256(hash_data).digest()[:32]

    # Compute u = y - (hash_value * r)
    hash_r = ensure_bytes(cp.scalar_multiply_scalar(hash_value, r))
    u = ensure_bytes(cp.scalar_subtraction(y, hash_r))

    # Compute v = z - (hash_value * s)
    hash_s = ensure_bytes(cp.scalar_multiply_scalar(hash_value, s))
    v = ensure_bytes(cp.scalar_subtraction(z, hash_s))

    return (hash_value, u, v, f_values)

def zkpok_verify(h_values: list[bytes], w_values: list[bytes], proof: tuple[bytes, bytes, bytes, list[bytes]]) -> bool:
    hash_value, u, v, f_values = proof

    # Compute M_i values for each index
    M_values = []
    for i in range(len(h_values)):
        # f_i^hash_value = point_multiply(hash_value, f_values[i])
        f_i_hash = ensure_bytes(cp.point_multiply(hash_value, f_values[i]))
        # h_i^u = point_multiply(u, h_values[i])
        h_i_u = ensure_bytes(cp.point_multiply(u, h_values[i]))
        # w_i^v = point_multiply(v, w_values[i])
        w_i_v = ensure_bytes(cp.point_multiply(v, w_values[i]))
        # med_i = f_i^hash_value + h_i^u
        med_i = ensure_bytes(cp.point_addition(f_i_hash, h_i_u))
        # M_i = med_i + w_i^v
        M_i = ensure_bytes(cp.point_addition(med_i, w_i_v))
        M_values.append(M_i)

    # Create hash input by concatenating f_values and M_values as bytes
    hash_data = b"".join(f_values + M_values)
    # Compute SHA-256 hash and truncate to 32 bytes
    hash_value_prime = hashlib.sha256(hash_data).digest()[:32]

    # Verification: check if recomputed hash matches the provided hash_value
    return hash_value_prime == hash_value

# ---------------------------------------------------------------------------
# Dealer commitment  cm = H1^r · H2^s  (same for all players)
# ---------------------------------------------------------------------------
def compute_cm(r: bytes, s: bytes) -> bytes:
    h1r = ensure_bytes(cp.point_multiply(r, H1_POINT))
    h2s = ensure_bytes(cp.point_multiply(s, H2_POINT))
    return ensure_bytes(cp.point_addition(h1r, h2s))


# ---------------------------------------------------------------------------
# Fiat–Shamir ZKPoK proving BOTH equations with the same (r,s)
# ---------------------------------------------------------------------------
def zkpok_proof_full(h_vals, w_vals, f_vals, cm, r, s):
    """
    Proves knowledge of r,s such that
        cm      = H1^r · H2^s             and
        f_i     = h_i^r · w_i^s  for all i
    """
    k_r, k_s = random_scalar(), random_scalar()          # nonce pair

    # Commitments
    A_cm = ensure_bytes(cp.point_addition(
              ensure_bytes(cp.point_multiply(k_r, H1_POINT)),
              ensure_bytes(cp.point_multiply(k_s, H2_POINT))))

    A_fis = [ensure_bytes(cp.point_addition(
                 ensure_bytes(cp.point_multiply(k_r, h_i)),
                 ensure_bytes(cp.point_multiply(k_s, w_i))))
             for h_i, w_i in zip(h_vals, w_vals)]

    # Fiat–Shamir challenge
    c_input = [cm] + f_vals + [A_cm] + A_fis
    c = hashlib.sha256(b"".join(c_input)).digest()[:32]

    # Responses
    z_r = ensure_bytes(cp.scalar_subtraction(
              k_r, ensure_bytes(ensure_bytes(cp.scalar_multiply_scalar(c, r)))))
    z_s = ensure_bytes(cp.scalar_subtraction(
              k_s, ensure_bytes(ensure_bytes(cp.scalar_multiply_scalar(c, s)))))

    return (c, z_r, z_s)        # cm and f_vals are public; no need to embed


def zkpok_verify_full(h_vals, w_vals, f_vals, cm, proof):
    c, z_r, z_s = proof

    # M_cm  = cm^c · H1^z_r · H2^z_s
    M_cm = ensure_bytes(cp.point_addition(
              ensure_bytes(cp.point_multiply(c, cm)),
              ensure_bytes(cp.point_addition(
                  ensure_bytes(cp.point_multiply(z_r, H1_POINT)),
                  ensure_bytes(cp.point_multiply(z_s, H2_POINT))))))

    # M_i   = f_i^c · h_i^z_r · w_i^z_s
    M_fis = [ensure_bytes(cp.point_addition(
                 ensure_bytes(cp.point_multiply(c, f_i)),
                 ensure_bytes(cp.point_addition(
                     ensure_bytes(cp.point_multiply(z_r, h_i)),
                     ensure_bytes(cp.point_multiply(z_s, w_i))))))
             for h_i, w_i, f_i in zip(h_vals, w_vals, f_vals)]

    c_prime = hashlib.sha256(b"".join([cm] + f_vals + [M_cm] + M_fis)).digest()[:32]
    return c_prime == c

def test_protocol_no_mask():
    print(f"==========TEST PROTOCOL NO MASK==========")
    print(f"Generated {n} shares with threshold {k}.")
    # Generate glist with t random points
    glist = gen_glist(k-1)
    
    # Generate x_values as distinct 32-byte little-endian integers (1 to n)
    x_values = [int.to_bytes(i + 1, 32, 'little') for i in range(n)]

    # Initialize hlist and wlist
    hlist = []
    wlist = [GENERATOR_POINT for _ in range(n)]
    
    process = psutil.Process(os.getpid())
    start = time.perf_counter()
    for i in range(n):
        h_i = compute_hp_no_mask(x_values[i], glist)
        hlist.append(h_i)
    stop = time.perf_counter()
    print(f"===> [1] Time to generate {n} hp = {(stop - start)/n:.3f}")

    # Generate random scalars r and s
    r = random_scalar()
    s = random_scalar()
    # Compute fplist = [fp_1, fp_2, ..., fp_n]
    start = time.perf_counter()
    fplist = compute_fp(hlist, wlist, r, s)
    stop = time.perf_counter()
    print(f"===> [2] Time to generate {n} fp = {stop - start:.3f}")
    
    # Generate proof using zkpok_proof
    start = time.perf_counter()
    proof = zkpok_proof(hlist, wlist, fplist, r, s)
    stop = time.perf_counter()
    print(f"===> [3] Time to generate a PV proof: {stop - start:.3f} seconds")
    
    # Verify proof using zkpok_verify
    start = time.perf_counter()
    is_valid = zkpok_verify(hlist, wlist, proof)
    stop = time.perf_counter()
    print(f"===> [4] Time to verify a PV proof: {stop - start:.3f} seconds - <Proof valid: {is_valid}>")
    """
    rec
    """
    # Compute shares for each player
    shares = []
    for i in range(n):
        shares.append((x_values[i], fplist[i]))


    # Reconstruct
    original_secret = ensure_bytes(cp.scalar_multiply(s))
    selected_shares = shares[:k]
    start = time.perf_counter()
    reconstructed_secret = reconstruct(selected_shares)
    stop = time.perf_counter()
    print(f"===> [7] Time to reconstruction = {stop-start:.3f}")
    if reconstructed_secret == original_secret:
        print("Test passed: Reconstructed secret matches g_0^s.")
    else:
        print("Test failed: Reconstructed secret does not match g_0^s.")
        print(f"Original secret (g_0^s): {original_secret.hex()}")
        print(f"Reconstructed secret: {reconstructed_secret.hex()}")

# ---------------------------------------------------------------------------
# ─── UPDATED TEST PROTOCOL (NO MASK VARIANT) ───────────────────────────────
# ---------------------------------------------------------------------------
def test_protocol_no_mask2():
    """
    • Builds h_i = f(x_i)  and  w_i = g
    • Dealer picks global r,s
    • Publishes
          cm  = H1^r · H2^s         (same for all players)
          ρ(x_i)=h_i^r · g^s        (called fp_i in code)
    • Produces / verifies the Fiat–Shamir proof that both
      equations use the *same* (r,s)
    • Reconstructs g^s from the first k shares
    """
    print("========== TEST PROTOCOL NO MASK ==========")
    print(f"Generated {n} shares with threshold {k}.")

    # 1) public polynomial points h_i  (= f(x_i))
    glist  = gen_glist(k - 1)
    x_vals = [int.to_bytes(i + 1, 32, 'little') for i in range(n)]
    hlist  = [compute_hp_no_mask(x, glist) for x in x_vals]

    # 2) single public base g for all w_i
    wlist = [GENERATOR_POINT] * n

    # 3) dealer’s randomness and constant commitment
    r, s = random_scalar(), random_scalar()
    cm   = compute_cm(r, s)

    # 4) ρ(x_i)  (called fp_i in earlier code)
    fplist = compute_fp(hlist, wlist, r, s)

    # 5) Fiat–Shamir ZK proof for both equations with SAME (r,s)
    proof = zkpok_proof_full(hlist, wlist, fplist, cm, r, s)
    proof_ok = zkpok_verify_full(hlist, wlist, fplist, cm, proof)
    print(f"   → Proof verifies? {proof_ok}")

    # 6) threshold reconstruction of  g^s
    shares  = list(zip(x_vals, fplist))     # (x_i , ρ(x_i))
    secret  = ensure_bytes(cp.scalar_multiply(s))  # g^s
    rec     = reconstruct(shares[:k])              # first k shares
    print(f"   → Reconstruction OK? {rec == secret}\n")


# Main execution
if __name__ == "__main__":
    print("Running Curve25519 property tests...")
    print("-" * 40)
    test_protocol_no_mask2()
    print("-" * 40)
    print("All tests completed!")
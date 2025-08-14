import secrets
import hashlib
import curve25519_python as cp
import time
import psutil
import os
IDENTITY_POINT = b'\x01' + b'\x00' * 31

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

def random_scalar():
    """Generate a random 32-byte scalar."""
    return secrets.token_bytes(32)

def random_point():
    """Generate a random point by multiplying the base point by a random scalar."""
    scalar = random_scalar()
    return ensure_bytes(cp.scalar_multiply(scalar))

# Core functions with ensured bytes
def gen_glist(t: int) -> list[bytes]:
    """
    Generate a list of t random points on Curve25519.
    Args:
        t: The number of random points to generate.
    Returns:
        A list of t 32-byte compressed EdwardsPoints as bytes.
    """
    return [ensure_bytes(cp.scalar_multiply(random_scalar())) for _ in range(t)]

def compute_hp(x: bytes, m: bytes, glist: list[bytes]) -> bytes:
    """
    Compute hp = (sum_{k=1}^t g_k^{x^k})^m where g_k are points in glist.
    Args:
        x: 32-byte scalar as bytes.
        m: 32-byte scalar as bytes.
        glist: List of 32-byte points (compressed EdwardsPoints) as bytes.
    Returns:
        32-byte compressed EdwardsPoint as bytes representing hp.
    """
    IDENTITY_POINT = b'\x01' + b'\x00' * 31
    h = IDENTITY_POINT
    x_pow = x
    for g in glist:
        p_k = ensure_bytes(cp.point_multiply(x_pow, g))
        h = ensure_bytes(cp.point_addition(h, p_k))
        x_pow = ensure_bytes(cp.scalar_multiply_scalar(x_pow, x))
    # print(f"==> There are {i} muls in compute_hp\n")
    hp = ensure_bytes(cp.point_multiply(m, h))
    return hp

def compute_hp_no_mask(x: bytes, glist: list[bytes]) -> bytes:
    hp = IDENTITY_POINT
    x_pow = x
    i = 0
    for g in glist:
        i = i+1
        p_k = ensure_bytes(cp.point_multiply(x_pow, g))
        hp = ensure_bytes(cp.point_addition(hp, p_k))
        x_pow = ensure_bytes(cp.scalar_multiply_scalar(x_pow, x))
    return hp

def compute_w(m: bytes):
    return ensure_bytes(cp.scalar_multiply(m))

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
    """
    Reconstruct the secret point from k shares using Lagrange interpolation.

    Args:
        shares: List of (x_j, y_j) tuples, where x_j is a scalar and y_j is a point.

    Returns:
        The reconstructed secret point as a 32-byte compressed Edwards point.
    """
    k = len(shares)
    IDENTITY_POINT = b'\x01' + b'\x00' * 31
    total_yj_delta = IDENTITY_POINT
    for j in range(k):
        t = 0
        xj_bytes, yj_bytes = shares[j]
        delta_bytes = b'\x01' + b'\x00' * 31
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
    """Generate zero-knowledge proof π.

    Args:
        h_values: List of public points (32-byte each) as bytes.
        w_values: List of public points (32-byte each) as bytes.
        f_values: List of committed points (32-byte each) as bytes.
        r: Secret scalar (32-byte) as bytes.
        s: Secret scalar (32-byte) as bytes.

    Returns:
        Tuple containing (hash_value, u, v, f_values), all as bytes except f_values as list[bytes].
    """
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
    """Verify zero-knowledge proof π.

    Args:
        h_values: List of public points (32-byte each) as bytes.
        w_values: List of public points (32-byte each) as bytes.
        proof: Tuple of (hash_value, u, v, f_values) from zkpok_proof.

    Returns:
        Boolean indicating whether the proof is valid.
    """
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

def zkpok_verify_one(h_value: bytes, w_value: bytes, proof: tuple[bytes, bytes, bytes, bytes]) -> bool:
    hash_value, u, v, f_value = proof

    # Compute intermediate points
    f_i_hash = ensure_bytes(cp.point_multiply(hash_value, f_value))  # f_i^hash_value
    h_i_u = ensure_bytes(cp.point_multiply(u, h_value))              # h_i^u
    w_i_v = ensure_bytes(cp.point_multiply(v, w_value))              # w_i^v

    # Compute M_i
    med_i = ensure_bytes(cp.point_addition(f_i_hash, h_i_u))         # med_i = f_i^hash_value + h_i^u
    M_i = ensure_bytes(cp.point_addition(med_i, w_i_v))              # M_i = med_i + w_i^v

    # Correctly concatenate f_value and M_i as bytes
    hash_data = f_value + M_i  # Direct concatenation of two bytes objects

    # Compute SHA-256 hash and truncate to 32 bytes
    hash_value_prime = hashlib.sha256(hash_data).digest()[:32]

    # Verification: check if recomputed hash matches the provided hash_value
    return hash_value_prime == hash_value
def zkpok_proof_one(h_value: bytes, w_value: bytes, f_value: bytes, r: bytes, s: bytes) -> tuple[bytes, bytes, bytes, bytes]:
    # Sample random scalars y and z as 32-byte objects
    y = secrets.token_bytes(32)
    z = secrets.token_bytes(32)

    # Compute commitments
    h_i_y = ensure_bytes(cp.point_multiply(y, h_value))
    w_i_z = ensure_bytes(cp.point_multiply(z, w_value))
    commitment = ensure_bytes(cp.point_addition(h_i_y, w_i_z))

    # Correctly concatenate f_value and commitment as bytes
    hash_data = f_value + commitment  # Direct concatenation of two bytes objects

    # Compute the hash value (SHA-256, truncated to 32 bytes)
    hash_value = hashlib.sha256(hash_data).digest()[:32]

    # Compute u = y - (hash_value * r)
    hash_r = ensure_bytes(cp.scalar_multiply_scalar(hash_value, r))
    u = ensure_bytes(cp.scalar_subtraction(y, hash_r))

    # Compute v = z - (hash_value * s)
    hash_s = ensure_bytes(cp.scalar_multiply_scalar(hash_value, s))
    v = ensure_bytes(cp.scalar_subtraction(z, hash_s))

    return (hash_value, u, v, f_value)
####################### TEST SPACE #######################

# Test functions
def test_scalar_invert():
    """
    Check if s1 * s1^{-1} == e, where e is the multiplicative identity (1).
    Note: The query mentions 's1 + s1^-1 == e', but in the scalar field, we check multiplication.
    """
    s1 = random_scalar()
    try:
        s1_inv = ensure_bytes(cp.scalar_inverse(s1))
    except ValueError as e:
        if "Cannot invert zero" in str(e):
            print("Scalar invert test skipped: s1 is zero")
            return
        else:
            raise
    identity = b'\x01' + b'\x00' * 31
    result = ensure_bytes(cp.scalar_multiply_scalar(s1, s1_inv))
    if result == identity:
        print("Scalar invert test passed: s1 * s1^{-1} == 1")
    else:
        print("Scalar invert test failed: s1 * s1^{-1} != 1")
def test_point_addition():
    """Check if p1 = p2 + p3 and p2 = p1 - p3."""
    p2 = random_point()
    p3 = random_point()
    p1 = ensure_bytes(cp.point_addition(p2, p3))
    p2_prime = ensure_bytes(cp.point_subtraction(p1, p3))
    if p2_prime == p2:
        print("Point addition test passed: p2 + p3 - p3 == p2")
    else:
        print("Point addition test failed: p2 + p3 - p3 != p2")
def test_point_multiplication():
    """Check if p1 = p2^{s1} and p2 = p1^{s1^{-1}}."""
    s1 = random_scalar()
    p2 = random_point()
    try:
        s1_inv = ensure_bytes(cp.scalar_inverse(s1))
    except ValueError as e:
        if "Cannot invert zero" in str(e):
            print("Point multiplication test skipped: s1 is zero")
            return
        else:
            raise
    p1 = ensure_bytes(cp.point_multiply(s1, p2))
    p2_prime = ensure_bytes(cp.point_multiply(s1_inv, p1))
    if p2_prime == p2:
        print("Point multiplication test passed: p2^{s1}^{s1^{-1}} == p2")
    else:
        print("Point multiplication test failed: p2^{s1}^{s1^{-1}} != p2")
def generate_shares(secret_scalar: bytes, n: int, k: int, t: int, r: bytes, s: bytes, glist: list[bytes]) -> list[tuple[bytes, bytes]]:
    x_values = [int.to_bytes(i + 1, 32, 'little') for i in range(n)]
    shares = []
    for x in x_values:
        start = time.perf_counter()
        m_i = random_scalar()
        hp_i = compute_hp(x, m_i, glist)
        w_i = compute_w(m_i)
        fp_i = compute_fp([hp_i], [w_i], r, s)[0]
        m_i_inv = ensure_bytes(cp.scalar_inverse(m_i))
        f_i = ensure_bytes(cp.point_multiply(m_i_inv, fp_i))
        shares.append((x, f_i))
        stop = time.perf_counter()
    print(f"generate 1 share = {stop-start}")
    return shares
def test_modified_reconstruction():
    s = random_scalar()
    original_secret = ensure_bytes(cp.scalar_multiply(s))
    print(f"Original secret (g_0^s): {original_secret.hex()}")
    # n, k = 256, 171
    t = k-1
    r = random_scalar()
    glist = gen_glist(t)
    shares = generate_shares(s, n, k, t, r, s, glist)
    print(f"Generated {n} shares with threshold {k}.")
    selected_shares = shares[:k]
    start = time.perf_counter()
    reconstructed_secret = reconstruct(selected_shares)
    stop = time.perf_counter()
    print(f"reconstruction = {stop-start}")
    print(f"Reconstructed secret: {reconstructed_secret.hex()}")
    if reconstructed_secret == original_secret:
        print("Test passed: Reconstructed secret matches g_0^s.")
    else:
        print("Test failed: Reconstructed secret does not match g_0^s.")
        print(f"Original secret (g_0^s): {original_secret.hex()}")
        print(f"Reconstructed secret: {reconstructed_secret.hex()}")
def test_zkp():
    # Example inputs (for illustration; actual values should be computed appropriately)
    glist = gen_glist(k-1)
    x_values = [int.to_bytes(i + 1, 32, 'little') for i in range(n)]
    m_values = [random_scalar() for i in range(n)]
    m_inv = []
    shares = []
    hlist = []
    wlist = []
    for i in range(n):
        hlist[i] = compute_hp(x_values[i], m_values[i], glist)
        wlist[i] = compute_w(m_values[i])
        m_inv[i] = ensure_bytes(cp.scalar_inverse(m_values[i]))

    r = secrets.token_bytes(32)
    s = secrets.token_bytes(32)
    # flist = compute_fp(hlist, wlist, r, s)
    fplist = compute_fp(hlist, wlist, r, s)



    # hlist = [ensure_bytes(cp.scalar_multiply(secrets.token_bytes(32))) for _ in range(n)]
    # wlist = [ensure_bytes(cp.scalar_multiply(secrets.token_bytes(32))) for _ in range(n)]
    # f_values = [ensure_bytes(cp.scalar_multiply(secrets.token_bytes(32))) for _ in range(122)]



    # Generate proof
    start = time.perf_counter()
    proof = zkpok_proof(hlist, wlist, flist, r, s)
    stop = time.perf_counter()
    print(f"===> Time to generate a proof: {stop-start}")


    # Verify proof
    start = time.perf_counter()
    is_valid = zkpok_verify(hlist, wlist, proof)
    stop = time.perf_counter()
    print(f"===> Time to verify a proof: {stop-start} - <Proof valid: {is_valid}>")
    
    for i in range(n):
        f_i = ensure_bytes(cp.point_multiply(m_inv[i], fplist[i]))
        shares.append((x_values[i], f_i))
        
def test_protocol():
    """
    Test the ZKP generation and verification with n players and t terms for compute_hp.

    Args:
        n: Number of players (default: 30).
        t: Number of terms for compute_hp (default: 5).
    """
    print(f"Generated {n} shares with threshold {k}.")
    # Generate glist with t random points
    glist = gen_glist(k-1)
    
    # Generate x_values as distinct 32-byte little-endian integers (1 to n)
    x_values = [int.to_bytes(i + 1, 32, 'little') for i in range(n)]
    
    # Generate m_values and m_inv for each player
    m_values = [random_scalar() for _ in range(n)]
    m_inv = [ensure_bytes(cp.scalar_inverse(m)) for m in m_values]
    
    # Initialize hlist and wlist
    hlist = []
    wlist = []
    
    process = psutil.Process(os.getpid())
    # start = time.perf_counter()
    # h_i = compute_hp(x_values[0], m_values[0], glist)
    # w_i = compute_w(m_values[0])
    # stop = time.perf_counter()   
    # hlist.append(h_i)
    # wlist.append(w_i)
    # print(f"===> [1] Time to generate 1 hp and 1 wp = {(stop - start):.3f}")
    for i in range(n):
        start = time.perf_counter()
        h_i = compute_hp(x_values[i], m_values[i], glist)
        w_i = compute_w(m_values[i])

        """
        n >= 512
        """
        # h_i = hlist[0]
        # w_i = wlist[0]
        ##
        hlist.append(h_i)
        wlist.append(w_i)
    stop = time.perf_counter()
    print(f"===> [1] Time to generate 1 hp and 1 wp = {(stop - start)/n:.3f}")

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
    mem_info1 = process.memory_info()
    print(f"RSS (Resident Set Size BEFORE): {mem_info1.rss / 1024 / 1024:.2f} MB")
    start = time.perf_counter()
    is_valid = zkpok_verify(hlist, wlist, proof)
    stop = time.perf_counter()
    mem_info2 = process.memory_info()
    print(f"RSS (Resident Set Size AFTER): {mem_info2.rss / 1024 / 1024:.2f} MB")
    print(f"===> [4] Time to verify a PV proof: {stop - start:.3f} seconds - <Proof valid: {is_valid}>")

    """
    test non pv
    """
    # Generate proof using zkpok_proof
    
    mem_info1 = process.memory_info()
    print(f"RSS (Resident Set Size BEFORE): {mem_info1.rss / 1024 / 1024:.2f} MB")
    start = time.perf_counter()
    for i in range(n):
        proof_one = zkpok_proof_one(hlist[i], wlist[i], fplist[i], r, s)
    stop = time.perf_counter()
    mem_info2 = process.memory_info()
    print(f"RSS (Resident Set Size AFTER): {mem_info2.rss / 1024 / 1024:.2f} MB")
    print(f"===> [5] Time to generate {n} proof (Non Public verifiability): {stop - start:.3f} seconds")

    # Verify proof using zkpok_verify
    start = time.perf_counter()
    is_valid_one = zkpok_verify_one(hlist[n-1], wlist[n-1], proof_one)
    stop = time.perf_counter()
    print(f"===> [6] Time to verify a non-PV proof: {stop - start:.3f} seconds - <Proof valid: {is_valid_one}>")
    """
    rec
    """
    
    # Compute shares for each player
    shares = []
    for i in range(n):
        # f_i = fp_i^{m_i^{-1}}
        f_i = ensure_bytes(cp.point_multiply(m_inv[i], fplist[i]))
        shares.append((x_values[i], f_i))


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

def test_protocol_no_mask():
    print(f"Generated {n} shares with threshold {k}.")
    # Generate glist with t random points
    glist = gen_glist(k-1)
    
    # Generate x_values as distinct 32-byte little-endian integers (1 to n)
    x_values = [int.to_bytes(i + 1, 32, 'little') for i in range(n)]

    # Initialize hlist and wlist
    hlist = []
    wlist = [IDENTITY_POINT for _ in range(n)]
    
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
    test non pv
    """
    # Generate proof using zkpok_proof
    
    mem_info1 = process.memory_info()
    print(f"RSS (Resident Set Size BEFORE): {mem_info1.rss / 1024 / 1024:.2f} MB")
    start = time.perf_counter()
    for i in range(n):
        proof_one = zkpok_proof_one(hlist[i], wlist[i], fplist[i], r, s)
    stop = time.perf_counter()
    mem_info2 = process.memory_info()
    print(f"RSS (Resident Set Size AFTER): {mem_info2.rss / 1024 / 1024:.2f} MB")
    print(f"===> [5] Time to generate {n} proof (Non Public verifiability): {stop - start:.3f} seconds")

    # Verify proof using zkpok_verify
    start = time.perf_counter()
    is_valid_one = zkpok_verify_one(hlist[n-1], wlist[n-1], proof_one)
    stop = time.perf_counter()
    print(f"===> [6] Time to verify a non-PV proof: {stop - start:.3f} seconds - <Proof valid: {is_valid_one}>")
    """
    rec
    """
    
    # Compute shares for each player
    shares = []
    for i in range(n):
        # f_i = fp_i^{m_i^{-1}}
        # f_i = ensure_bytes(cp.point_multiply(m_inv[i], fplist[i]))
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


# Main execution
if __name__ == "__main__":
    print("Running Curve25519 property tests...")
    print("-" * 40)
    # test_scalar_invert()
    # test_point_addition()
    # test_point_multiplication()
    # test_reconstruction()
    test_modified_reconstruction()
    # test_protocol_no_mask()
    test_protocol()
    print("-" * 40)
    print("All tests completed!")
import secrets, time, gc, tracemalloc

# Start heap tracing
tracemalloc.start()

def heap_kib():
    cur, peak = tracemalloc.get_traced_memory()
    return cur // 1024, peak // 1024

def banner():
    gc.collect()
    cur, peak = heap_kib()
    print(f"   [heap] {cur:6d} KiB  (peak {peak} KiB)")

# Normalize to bytes
def E(v):
    return v if isinstance(v, (bytes, bytearray)) else bytes(v)

# Global counters
CNT = RB_CNT = 0
RB_TIME = 0.0

def bump():
    global CNT
    CNT += 1

def rb_bump(d):
    global RB_CNT
    RB_CNT += d

def reset():
    global CNT
    CNT = 0

def rb_reset():
    global RB_CNT, RB_TIME
    RB_CNT = 0
    RB_TIME = 0.0

def rand():
    return secrets.token_bytes(32)
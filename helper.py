def compute_cost_Trace(n, t, f):
    return (n - f) * ((t + 5) * (t - f) + t) + f * (t + 5) * (t - f)

def compute_cost_TrVer(n, t, f):
    return 2 * n * t  + 3 + f * ((t + 5) * (t - f))

# def compute_cost_TrVer(n, t, f):
#     return 2 * n * (t+1) + f * (t + 5)

    # 2n(t+1) + f(t+5)

if __name__ == "__main__":
    for n, t, f in [(32, 17, 11), (64, 33, 22), (128, 65, 43), (256, 129, 86)]:
        print(f"compute_cost_Trace({n},{t},{f}) = {compute_cost_Trace(n, t, f)}")
        print(f"compute_cost_TrVer({n},{t},{f}) = {compute_cost_TrVer(n, t, f)}")
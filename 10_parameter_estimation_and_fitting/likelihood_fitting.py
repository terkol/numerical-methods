import numpy as np
import matplotlib.pyplot as plt
import os

def L(data, x1, x2):
    return np.prod([(1/(np.sqrt(2*np.pi)*i[2])*np.exp(-((i[1]-(x1*i[0]+x2))**2)/(2*i[2]*i[2]))) for i in data]) # Likelihood function for given x1 and x2

def metropolis(x0, data, steps=100000, burnin=10000, report=True):
    # Markov Chain Monte Carlo sampler (intentionally uses raw likelihood product)
    x1, x2 = x0

    # Proposal scales (keep your original idea; optionally use abs(...) to avoid freezing)
    max_d1 = x1 * 0.05
    max_d2 = x2 * 0.05

    rands1 = np.random.rand(steps) * 2 - 1
    rands2 = np.random.rand(steps) * 2 - 1
    u = np.random.rand(steps)

    # Instrumentation counters
    stats = {
        "L_zero": 0,
        "L_inf": 0,
        "L_nan": 0,
        "q_nan": 0,
        "accepted": 0,
    }

    lh = L(data, x1, x2)
    if lh == 0.0:
        stats["L_zero"] += 1
    elif np.isinf(lh):
        stats["L_inf"] += 1
    elif np.isnan(lh):
        stats["L_nan"] += 1

    chain = np.empty((steps, 2), dtype=float)

    for i in range(steps):
        new_x1 = x1 + rands1[i] * max_d1
        new_x2 = x2 + rands2[i] * max_d2

        new_lh = L(data, new_x1, new_x2)

        # Track numerical pathologies in likelihood itself
        if new_lh == 0.0:
            stats["L_zero"] += 1
        elif np.isinf(new_lh):
            stats["L_inf"] += 1
        elif np.isnan(new_lh):
            stats["L_nan"] += 1

        # Compute acceptance ratio q = new_lh / lh, but handle 0/0 and NaNs explicitly
        if lh == 0.0:
            # If lh underflowed to 0, the ratio is undefined; this is part of the instability.
            # Policy: if new_lh > 0, treat as "inf"; else treat as NaN -> reject.
            if new_lh > 0.0 and np.isfinite(new_lh):
                q = np.inf
            else:
                q = np.nan
        else:
            q = new_lh / lh

        if not np.isfinite(q):
            # Count NaN/inf ratios; reject NaNs, accept infs
            if np.isnan(q):
                stats["q_nan"] += 1
                accept = False
            else:
                accept = True  # q = +inf
        else:
            accept = (q >= 1.0) or (u[i] < q)

        if accept:
            x1, x2, lh = new_x1, new_x2, new_lh
            stats["accepted"] += 1

        chain[i, 0] = x1
        chain[i, 1] = x2

    if report:
        acc_rate = stats["accepted"] / steps
        print("=== Raw-likelihood instability counters ===")
        print(f"steps: {steps}, burn-in: {burnin}, acceptance rate: {acc_rate:.3f}")
        print(f"L == 0  count: {stats['L_zero']}")
        print(f"L == inf count: {stats['L_inf']}")
        print(f"L == NaN count: {stats['L_nan']}")
        print(f"q == NaN count: {stats['q_nan']}")
        print("==========================================")

    return chain[burnin:]

if __name__ == "__main__":
    data = np.loadtxt(os.path.dirname(os.path.dirname(__file__))+"\\data\\raw\\linear_data_errors.txt")
    print(data)
    x0 = (2,5)

    mc = metropolis(x0, data, report=True)
    x1 = mc[:,0]
    x2 = mc[:,1]

    plt.figure(figsize=(12,6))
    plt.subplot(1,2,1)
    plt.hist(x1, 100, (1.5,2.5), density=True)
    plt.title("x1")
    plt.grid()
    plt.subplot(1,2,2)
    plt.hist(x2, 100, (4,8), density=True)
    plt.title("x2")
    plt.grid()
    plt.show()

    plt.errorbar(data[:,0],data[:,1], yerr=data[:,2], fmt="o", capsize=5)
    plt.xlabel("t")
    plt.ylabel("y")
    plt.grid()
    plt.show()


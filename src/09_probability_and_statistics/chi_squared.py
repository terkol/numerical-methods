import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import chi2
import time

def gen_points(num, num_bins):  # Generate points using np.random.randint and bin them
    bins = np.bincount(rands, minlength=num_bins)
    rands = np.random.randint(0,num_bins,num)   # Generate random numbers between 0 and num_bins
    for r in rands:     # For each random number
        bins[r] += 1    # Increment bin value for matching index
    return bins         # Return array of bins

def chi_squared(bins):
    return sum([((np.mean(bins)-i)**2)/np.mean(bins) for i in bins]) # Return chi squared from given bins

def bad_lcg(x, a=1103515245, c=12345, m=65536): # Bad lc generator
    while True:
        x = (a*x+c)%m
        yield x
        
def gen_bad_points(gen, num, num_bins, m=65536): # Generate points using bad_lcg
    bins = np.bincount(rands, minlength=num_bins)
    rands = [int((next(gen)/m)*num_bins) for _ in range(num)]
    for r in rands:
        bins[r] += 1
    return bins

if __name__ == "__main__":
    plt.figure(figsize=(8,6))
    gen = bad_lcg(1) # Initialize bad lc generator

    start = time.time()
    chi_scores = [chi_squared(gen_bad_points(gen, 10000,10)) for _ in range(50000)]
    stop = time.time()
    print(stop-start)

    x = np.linspace(0, 40, 100)
    plt.plot(x, chi2.pdf(x, df=9))
    plt.hist(chi_scores,100,(0,40), density=True)
    plt.title("Using bad LCG")
    plt.xlabel("χ2")
    plt.ylabel("Probability")
    plt.grid()

    plt.figure(figsize=(8,6))
    start = time.time()
    chi_scores = [chi_squared(gen_points(gen, 10000,10)) for _ in range(50000)]
    stop = time.time()
    print(stop-start)

    x = np.linspace(0, 40, 100)
    plt.plot(x, chi2.pdf(x, df=9))
    plt.hist(chi_scores,100,(0,40), density=True)
    plt.title("Using a good RNG")
    plt.xlabel("χ2")
    plt.ylabel("Probability")
    plt.grid()

    plt.show()
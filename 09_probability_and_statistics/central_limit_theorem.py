import numpy as np
import matplotlib.pyplot as plt

def a(K):   # Uniform distribution between -1/2 and 1/2
    return np.random.random(K)-0.5 

def b(K):   # Sample inverse CDF of Laplace distribution
    samples = np.random.random(K)
    for i,v in enumerate(samples):
        if v < 0.5:
            samples[i] = np.log(2*v)
        else:
            samples[i] = -1*np.log(2*(1-v))
    return samples

def c(K):   # Sample inverse CDF of the Lorentz distribution
    return np.tan(np.pi*np.random.random(K)-0.5)

if __name__ == "__main__":
    for K in [3,10,30,100]: # For each of given K values
        z = []
        for _ in range(10000):              
            samples = c(K)                                              # Choose function 
            z.append(np.mean(samples))      # Add mean of samples to z
        plt.hist(z, bins=100, density=True) # Plot histogram of means
        x = np.linspace(min(z), max(z), 200)
        mu, sigma = np.mean(z), np.std(z)
        norm = (1/(sigma*np.sqrt(2*np.pi)))*np.exp(-(x-mu)**2/(2*sigma**2)) # Normal distribution
        plt.plot(x, norm, label="Normal fit")                    # Overlay normal distribution
        plt.title(f"K={K}")
        plt.legend()
        plt.show()

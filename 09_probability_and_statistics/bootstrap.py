import numpy as np
import matplotlib.pyplot as plt

if __name__ == "__main__":
    exps = 50
    B = 10000
    hits = 0
    for _ in range(exps):
        obs = np.random.exponential(5, 100) 
        means = []                 
        for b in range(B):
            resample = np.random.choice(obs, size=100, replace=True)  
            means.append(np.mean(resample))                      
        ci_low, ci_high = np.percentile(means, [2.5, 97.5])       
        if ci_low < 5 and ci_high > 5: 
            hits += 1
    print(f"{hits}/{exps} confidence intervals contained the true mean")
    plt.hist(means, bins=100, density=True)
    plt.title("Bootstrap distribution from the last experiment")
    plt.show()
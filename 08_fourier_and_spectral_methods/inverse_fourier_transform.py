import numpy as np
import matplotlib.pyplot as plt

if __name__ == "__main__":
    N = 128
    d = 2*np.pi/N
    u0 = []
    for j in range(N):
        u0.append(np.abs(np.cos(d*j)))
    u0_hat = np.fft.fft(u0)

    k = np.fft.fftfreq(N, d=d) * 2*np.pi

    times = [0, 0.2, 0.5, 10]
    solutions = {}
    for t in times:
        u_hat = u0_hat * np.exp(-(k**2) * t)
        solutions[t] = np.fft.ifft(u_hat).real

    x = np.linspace(0, 2*np.pi, N, endpoint=False)
    for t in times:
        plt.plot(x, solutions[t], label=f"t = {t}")
        
    plt.xlabel("x")
    plt.ylabel("u(x,t)")
    plt.legend()
    plt.grid()
    plt.show()
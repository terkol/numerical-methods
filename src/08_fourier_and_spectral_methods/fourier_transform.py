import numpy as np
import matplotlib.pyplot as plt
import pyfftw
from scipy.signal import spectrogram
from scipy.signal import find_peaks
import os

def fundamental_frequency():     # Function for finding the main frequency of the signal
    sig = signal-np.mean(signal) # Subtract mean
    N = len(sig)                 # Length of array
    freqs = np.fft.fftfreq(N, d=dt) # Find frequencies
    a = pyfftw.empty_aligned(N, dtype='complex128') # Initialize arrays for fftw
    b = pyfftw.empty_aligned(N, dtype='complex128')
    a.real = sig   
    fft_object = pyfftw.FFTW(a, b)  # Initialize fftw with parameters
    fft_object()    # Call fftw
    mag = np.abs(b)
    peak = np.argmax(mag[:N//2])          # restrict to nonnegative frequencies
    print("Peak Frequency:", freqs[:N//2][peak], "Hz")

def higher_harmonics(): # Work in progress
    sig = signal-np.mean(signal)
    N = len(sig)
    freqs = np.fft.fftfreq(N, d=dt)
    a = pyfftw.empty_aligned(N, dtype='complex128')
    a.real = sig
    b = pyfftw.empty_aligned(N, dtype='complex128')
    fft_object = pyfftw.FFTW(a, b)
    fft_object()
    fft_magnitude = np.abs(b[:N//2])
    freqs = np.fft.fftfreq(N, d=dt)
    peaks, properties = find_peaks(fft_magnitude, height=0.001*np.max(fft_magnitude), distance=1000)
    print("Peak indices:", peaks)
    print("Peak frequencies and heights:", np.array(list(zip(freqs[peaks], properties['peak_heights']))))

def plot_spectrogram(): # Plots spectrogram of signal
    f, t, Sxx = spectrogram(signal, fs=fs, nperseg=4096, window='hann') # High resolution spectrogram using Hanning window
    plt.pcolormesh(t, f, 10*np.log10(Sxx+1e-30))
    plt.ylim(0,4000) # Don't plot highest overtones as they aren't really visible anyway
    plt.ylabel("Hz")
    plt.xlabel("t")
    plt.show()

if __name__ == "__main__":
    data = np.loadtxt(os.path.dirname(os.path.dirname(os.path.dirname(__file__)))+'\\data\\sample1.dat', skiprows=1) # Load data
    time = data[:,0]    # Extract first column
    signal = data[:,1]  # Extract second column
    dt = time[1]-time[0]    # Find time interval
    fs = 1.0/dt         # Compute frequency of measurements

    fundamental_frequency()
    plot_spectrogram()


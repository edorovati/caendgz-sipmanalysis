# plotter.py
import matplotlib.pyplot as plt

class Plotter:
    def __init__(self):
        pass

    def plot_waveform(self, waveform, title, ylabel, label, color="blue"):
        plt.figure(figsize=(10, 5))
        plt.plot(waveform, label=label, color=color)
        plt.xlabel("Sample")
        plt.ylabel(ylabel)
        plt.title(title)
        plt.legend()
        plt.grid(True)
        plt.show()

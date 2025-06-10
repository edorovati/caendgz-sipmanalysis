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
        
    def plot_multiple_waveforms(self, waveforms, labels, title, xlabel, ylabel):
        for wf, label in zip(waveforms, labels):
            plt.plot(wf, label=label)
        plt.title(title)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.legend()
        plt.grid(True)
        plt.show()

    def plot_histogram(self, data, title, xlabel, ylabel, bins=50, color="blue"):
        plt.figure(figsize=(8, 6))
        plt.hist(data, bins=bins, color=color, alpha=0.7)
        plt.title(title)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.grid(True)
        plt.tight_layout()
        plt.show()

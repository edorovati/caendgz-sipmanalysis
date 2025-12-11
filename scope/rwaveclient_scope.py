import sys
import numpy as np
import threading
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from PyQt5.QtWidgets import (
    QApplication, QVBoxLayout, QWidget, QPushButton, QLabel, QCheckBox,
    QHBoxLayout, QSpacerItem, QSizePolicy
)
from PyQt5.QtCore import Qt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT
sys.path.append("/eu/caen-dt5742b/python/")
from rwave import rwaveclient

host = 'localhost'
port = 30001

class OscilloscopeApp(QWidget):
    def __init__(self):
        super().__init__()
        self.host = host
        self.port = port
        self.current_frequency = 750
        self.channel_visibility = {0: True, 1: True}
        self.latest_data = None
        self.current_frame = 0
        self.is_running = True
        self.correction_enabled = False
        self.rcorrection_enabled = False
        self.first_cell = 0 
        self.slopes = np.zeros(1024) 
        self.intercepts = np.zeros(1024) 
        self.initUI()
        self.init_plot()

    # ---------------------- Digitizer Control ----------------------
    def set_frequency(self, freq):
        self.current_frequency = freq
        threading.Thread(target=self.configure_digitizer, daemon=True).start()
    
    def configure_digitizer(self):
        with rwaveclient(self.host, self.port, verbose=True) as rwc:
            if rwc is not None:
                rwc.send_cmd(f'sampling {self.current_frequency}')
                rwc.send_cmd('grmask 0x1')
                rwc.send_cmd('chmask 0x0003')
                rwc.send_cmd('correction on' if self.correction_enabled else 'correction off')
    
    def toggle_channel(self, ch, state):
        self.channel_visibility[ch] = state == Qt.Checked
    
    def toggle_correction(self, state):
        self.correction_enabled = state == Qt.Checked
        threading.Thread(target=self.configure_digitizer, daemon=True).start()
    
    def toggle_rcorrection(self, state):
        self.rcorrection_enabled = state == Qt.Checked
        threading.Thread(target=self.update_plot_correction, daemon=True).start()

    # ---------------------- UI Setup ----------------------
    def initUI(self):
        # Main layout: canvas + side panel
        main_layout = QHBoxLayout()
        side_layout = QVBoxLayout()

        self.setWindowTitle("Digitizer Control")
        self.setGeometry(100, 100, 1000, 600)

        # Start Acquisition button
        self.button = QPushButton('Start Acquisition')
        self.button.setStyleSheet("background-color: #4CAF50; color: white; font-size: 12px; padding: 6px;")
        self.button.clicked.connect(self.startAcquisition)
        side_layout.addWidget(self.button)

        # Frequency buttons in one row
        freq_layout = QHBoxLayout()
        for freq in [750, 1000, 2500, 5000]:
            btn = QPushButton(f"{freq} MHz", self)
            btn.setStyleSheet("background-color: #008CBA; color: white; font-size: 12px; padding: 6px;")
            btn.clicked.connect(lambda _, f=freq: self.set_frequency(f))
            freq_layout.addWidget(btn)
        side_layout.addLayout(freq_layout)

        # Status label
        self.status_label = QLabel("Status: Ready")
        side_layout.addWidget(self.status_label)

        # Channel checkboxes
        for ch in [0, 1]:
            checkbox = QCheckBox(f'Channel {ch}', self)
            checkbox.setChecked(self.channel_visibility[ch])
            checkbox.stateChanged.connect(lambda state, ch=ch: self.toggle_channel(ch, state))
            side_layout.addWidget(checkbox)

        # Correction checkboxes
        self.correction_checkbox = QCheckBox('Enable Correction', self)
        self.correction_checkbox.stateChanged.connect(self.toggle_correction)
        side_layout.addWidget(self.correction_checkbox)

        self.rcorrection_checkbox = QCheckBox('Enable RCorrection', self)
        self.rcorrection_checkbox.stateChanged.connect(self.toggle_rcorrection)
        side_layout.addWidget(self.rcorrection_checkbox)

        # Navigation buttons
        self.next_frame_button = QPushButton('Next Frame', self)
        self.next_frame_button.clicked.connect(self.next_frame)
        side_layout.addWidget(self.next_frame_button)

        self.stop_frame_button = QPushButton('Stop Frame', self)
        self.stop_frame_button.clicked.connect(self.stop_frame)
        side_layout.addWidget(self.stop_frame_button)

        side_layout.addSpacerItem(QSpacerItem(20, 40, QSizePolicy.Minimum, QSizePolicy.Expanding))

        # Canvas + toolbar
        canvas_layout = QVBoxLayout()
        self.canvas = FigureCanvas(plt.figure(facecolor='black'))
        self.toolbar = NavigationToolbar2QT(self.canvas, self)
        canvas_layout.addWidget(self.toolbar)
        canvas_layout.addWidget(self.canvas)

        # Add layouts to main layout
        main_layout.addLayout(side_layout)
        main_layout.addLayout(canvas_layout, stretch=1)  # canvas occupies remaining space

        self.setLayout(main_layout)
        self.canvas.mpl_connect("motion_notify_event", self.toggle_animation)

    # ---------------------- Acquisition ----------------------
    def startAcquisition(self):
        self.button.setEnabled(False)
        self.status_label.setText("Status: Acquiring Data...")
        self.latest_data = self.acquireData()
        self.button.setEnabled(True)
        self.status_label.setText("Status: Ready")
        self.ani.event_source.start()
    
    def acquireData(self):
        with rwaveclient(self.host, self.port, verbose=True) as rwc:
            if rwc is None:
                return None
            rwc.send_cmd("start")
            rwc.send_cmd('swtrg 1024')
            rwc.send_cmd('readout')
            rwc.send_cmd('download')
            data = rwc.download()
            rwc.send_cmd('stop')
        return data

    # ---------------------- Plotting ----------------------
    def init_plot(self):
        self.ax = self.canvas.figure.add_subplot(111, facecolor='black')
        self.x_data = np.arange(1024)
        self.lines = {}
        for ch, color in zip([0, 1], ['yellow', 'cyan']):
            self.lines[ch], = self.ax.plot(self.x_data, np.zeros(1024), label=f'ch-{ch}', color=color)
        self.ax.set_ylim(0, 4096)
        self.ax.set_xlabel("Samples", color='white')
        self.ax.set_ylabel("Amplitude", color='white')
        self.ax.set_title("Oscilloscope Data", color='white')
        self.ax.grid(color='gray')
        self.ax.legend()
        self.ax.tick_params(axis='both', colors='white')
        self.canvas.figure.tight_layout()
        self.ani = animation.FuncAnimation(self.canvas.figure, self.update_plot, interval=100, blit=False, repeat=True)
        self.ani.event_source.stop()
    
    def toggle_animation(self, event):
        if self.toolbar.mode:
            self.ani.event_source.stop()
        else:
            self.ani.event_source.start()
    
    def next_frame(self):
        if self.latest_data is not None:
            self.current_frame = (self.current_frame + 1) % len(self.latest_data)
            self.update_plot(None)
    
    def stop_frame(self):
        self.ani.event_source.stop()
        self.is_running = False
    
    def update_plot(self, frame):
        if self.latest_data is not None:
            event = self.latest_data[self.current_frame]
            for ch in [0, 1]:
                if self.channel_visibility[ch]:
                    self.lines[ch].set_ydata(event[ch]["waveform"])
            self.ax.set_title(f'Frame {self.current_frame}', color='white')
            self.current_frame = (self.current_frame + 1) % len(self.latest_data)
            self.canvas.flush_events()

    # ---------------------- Resize handling ----------------------
    def resizeEvent(self, event):
        super().resizeEvent(event)

        # Canvas dynamic size: maintain 4:3, fit in available space
        available_width = self.canvas.parent().width()
        available_height = self.canvas.parent().height()

        new_width = min(available_width, int(available_height * 4 / 3))
        new_height = int(new_width * 3 / 4)

        self.canvas.resize(new_width, new_height)
        self.canvas.draw_idle()


if __name__ == '__main__':
    app = QApplication(sys.argv)
    window = OscilloscopeApp()
    window.show()
    sys.exit(app.exec_())

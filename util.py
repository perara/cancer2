import os
import warnings

from scipy.signal import find_peaks, peak_widths, find_peaks_cwt, savgol_filter
import re
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.signal._peak_finding_utils import PeakPropertyWarning
import pandas
sns.set()
warnings.simplefilter("ignore", PeakPropertyWarning)

HRR_PARAM_LOW = .60
HRR_PARAM_HIGH = .75
PEAK_WIDTH_LOW = 40  # 60 seconds
PEAK_WIDTH_HIGH = 250
TEST_MIN_POINTS = 10
END_HRR_STRIP = .40
DPI = 100
SAVE_PDF = False
HRR_LINES = [
    [.40, "#ff4000"],
    [.50, "#ff8000"],
    [.70, "#bfff00"],
    [.75, "#40ff00"],
    [.80, "#00ff80"],
    [.85, "#00bfff"],
    [.90, "#0040ff"],
    [.95, "#8000ff"],
    [1.0, "#ff00ff"],
]

HRR_FILL_COLOR = [
    [.90, 1.0, "#ff0000"],
    [.80, .90, "#dc4806"],
    [.50, .80, "#608d45"],
    [.40, .50, "#069adc"],
    [.0, .40, "#a7adba"]


]



class HRM:

    def __init__(self, path):
        self._path = path

        self._hr_data = None
        self._hr_interval = None
        self._max_hr = None
        self._rest_hr = None

        self._invalid = None
        self._invalids = []

        self.invalid = False
        self.rest_hr = None
        self.max_hr = None

        self._load()

    def validate(self):
        if len(self._hr_data) < TEST_MIN_POINTS:
            return False

        max_hr = self._hr_data.max()
        if max_hr < self.get_hrr_hr(HRR_PARAM_LOW):
            return False

        return True

    def _load_hr_data(self, _d):
        hr_data_index = _d.find('[HRData]')
        self._hr_data = np.array(list(map(float, re.findall(r'\d+', _d[hr_data_index:-1]))), dtype=np.int)
        strip_at = len(self._hr_data)

        try:
            peaks, _ = self.get_peaks_1()
            last_peak = peaks[-1]

            for i in range(last_peak, len(self._hr_data)):
                y = self._hr_data[i]
                if self.get_hrr(y) <= (END_HRR_STRIP*100):
                    strip_at = i
                    break
        except:
            pass

        self._hr_data = self._hr_data[:strip_at]

    def _load_hr_interval_value(self, _d):
        try:
            self._hr_interval = int(re.search(r'Interval=([0-9]+)', _d).group(1))
        except Exception as e:
            print("Replace this...", e)
            self._invalid = True
            self._invalids.append("Could not load hr_interval")

    def _load_hr_interval(self):
        """Load hr interval for series."""
        # TODO unscaled.
        self._hr_data_times = np.array([x for x in range(len(self._hr_data))], dtype=np.int)
        self._hr_data_times = self._hr_data_times

    def _load_max_hr(self, _d):
        self._max_hr = int(re.search(r'MaxHR=([0-9]+)', _d).group(1))

    def _load_rest_hr(self, _d):
        self._rest_hr = int(re.search(r'RestHR=([0-9]+)', _d).group(1))

    """def _parse_params(self):
        hr_data_start_idx = self.content.find('[Params]') + 8
        hr_data_end_idx = self.content.find('[Note]')
        hr_params_section_arr = self.content[hr_data_start_idx:hr_data_end_idx].strip().split("\n")
        hr_params_section = {
            v.split("=")[0]: v.split("=")[1] for v in hr_params_section_arr
        }

        self.rest_hr = float(hr_params_section["RestHR"])
        self.max_hr = float(hr_params_section["MaxHR"])"""

    def _load(self):

        with open(self._path, 'rb') as f:
            data = str(f.read())

        self._load_max_hr(data)
        self._load_rest_hr(data)
        self._load_hr_interval_value(data)
        self._load_hr_data(data)
        self._load_hr_interval()

    def get_hrr(self, hr):
        return (hr - self._rest_hr) / (self._max_hr - self._rest_hr) * 100.0

    def get_hrr_hr(self, hhr):
        return (hhr * (self._max_hr - self._rest_hr)) + self._rest_hr

    def get_peak_times(self, peaks):
        return np.multiply(peaks[1], self._hr_interval)

    def get_duration(self):
        """In seconds."""
        return len(self._hr_data) * self._hr_interval

    def smooth_filter(self):
        return savgol_filter(self._hr_data, 7, 4)

    def get_peaks_1(self):
        """peaks, _ = find_peaks(self.smooth_filter(),
                              height=self.get_hrr_hr(HRR_PARAM_LOW),
                              prominence=12,
                              width=[
                                  float(PEAK_WIDTH_LOW / self._hr_interval),
                                  float(PEAK_WIDTH_HIGH / self._hr_interval)
                              ]
                              )"""
        x = self.smooth_filter()
        import peakutils

        def find_bottom(data, peak, reverse=False):
            LOWEST = peak
            iter = range(peak, 0, -1) if reverse else range(peak, len(data))
            for i in iter:
                if data[i] > data[LOWEST]:
                    return LOWEST
                LOWEST = i
            return LOWEST

        peaks = np.array(peakutils.indexes(x))
        peaks_start = np.array([find_bottom(x, peak) for peak in peaks])
        peaks_end = np.array([find_bottom(x, peak, reverse=True) for peak in peaks])
        peaks_width = np.abs(np.array([x[i] for i in peaks_start]) - np.array([x[i] for i in peaks_end])) # TODO wrong?
        peaks_height = np.array([x[peaks[i]] - max(x[peaks_start[i]], x[peaks_end[i]]) for i in range(len(peaks))])

        df = pandas.DataFrame({
            "peaks": peaks,
            "peaks_width": peaks_width * self._hr_interval,
            "peaks_height": peaks_height
        })

        """Width constraint."""
        df = df.loc[(df['peaks_width'] >= PEAK_WIDTH_LOW) & (df['peaks_width'] <= PEAK_WIDTH_HIGH)]

        """Height constraint."""
        df = df.loc[(df['peaks_height'] >= 12)]



        _ = {}
        return df["peaks"].values, _

    def get_peaks(self):

        peaks, _ = self.get_peaks_1()

        # print(int(90 / self._hr_interval))
        hr_peaks = self._hr_data[peaks]
        peak_times = np.array([self._hr_data_times[x] for x in peaks])
        peak_w = peak_widths(self._hr_data, peaks, rel_height=0.5)

        return hr_peaks, peak_times, peak_w

    def get_values_over_hrr(self, hrr=0):
        bpm = self.get_hrr_hr(hrr)
        data = np.copy(self._hr_data)
        data[self._hr_data < bpm] = 0
        return data

    def plot(self, ax=None, save=False, save_location=None, save_name=None):
        if save is True and (save_location is None or save_name is None):
            raise RuntimeError("You must set a save location AND save_name")

        fig = None
        if ax is None:
            fig, ax = plt.subplots(1, 1)

        """Plot background."""
        for spec in HRR_FILL_COLOR:
            _low = self.get_hrr_hr(spec[0])
            _high = self.get_hrr_hr(spec[1])
            _color = spec[2]
            plt.fill_between(np.multiply(self._hr_data_times, self._hr_interval), _low, _high, color=_color, alpha=0.5)


        """Plot horizontal marker lines"""
        for _hrrs in HRR_LINES:
            plt.axhline(self.get_hrr_hr(_hrrs[0]), color=_hrrs[1], linewidth=1.0, label="HRR-%s" % _hrrs[0], linestyle="--")

        plt.axhline(self.get_hrr_hr(HRR_PARAM_LOW), color="C2", linewidth=3.0, label="LOW-%s" % HRR_PARAM_LOW, linestyle="--")
        plt.axhline(self.get_hrr_hr(HRR_PARAM_HIGH), color="C3", linewidth=3.0, label="HIGH-%s" % HRR_PARAM_HIGH, linestyle="--")



        """Plot lines."""
        xs = np.multiply(self._hr_data_times, self._hr_interval)
        ys = self.smooth_filter()
        sns.lineplot(xs, ys, color="red", linestyle="-", linewidth=1.0)
        #sns.lineplot(xs, self._hr_data, color="black")

        """Fill under line"""
        ax.fill_between(xs, ys, self.get_hrr_hr(.0), interpolate=True, color='#E0FFFF', alpha=0.9)

        if True:
            peaks, peak_times, peak_w = self.get_peaks()

            for p, pt in zip(peaks, peak_times):
                is_low = self.get_hrr(p) < (HRR_PARAM_HIGH*100)
                if is_low:
                    ax.plot(np.multiply(pt, self._hr_interval), p, "s", color="#ffa500")
                else:
                    ax.plot(np.multiply(pt, self._hr_interval), p, "s", color="#ff0000")

            """Scale width of hlines properly"""
            """Plot width of peak"""
            hl = peak_w[1:]
            a = hl[0]
            b = np.multiply(hl[1], self._hr_interval)
            c = np.multiply(hl[2], self._hr_interval)
            plt.hlines(a, b, c, color="C1")

            # results_half = peak_widths(self.rri, peaks_to_time, rel_height=0.5)

            # plt.hlines(*results_full[1:], color="C3")

        ax.set(xlabel='Time (s)', ylabel='Heart Rate (bps)')
        ax.margins(0.0)
        # plt.show(block=False)

        plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        plt.title(os.path.basename(self._path))
        plt.tight_layout()
        plt.savefig(os.path.join(save_location, save_name), bbox_inches="tight", pad_inches=0, dpi=DPI)
        if SAVE_PDF:
            plt.savefig(os.path.join(save_location, save_name) + ".pdf", bbox_inches="tight", pad_inches=0, dpi=DPI)
        plt.close(fig)

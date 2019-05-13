import os
import warnings
import matplotlib.patches as mpatches
from scipy.signal import find_peaks, peak_widths, find_peaks_cwt, savgol_filter, peak_prominences
import re
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.signal._peak_finding_utils import PeakPropertyWarning
import pandas
from scipy.interpolate import interp1d

sns.set()
warnings.simplefilter("ignore", PeakPropertyWarning)

METHOD = 1

HRR_PARAM_LOW = .60
HRR_PARAM_HIGH = .80
PEAK_WIDTH_LOW = 40  # 60 seconds
PEAK_WIDTH_HIGH = 500
TEST_MIN_POINTS = 22
PEAK_HEIGHT = 10
END_HRR_STRIP = .39
DPI = 100
SAVE_PDF = False
HRR_LINES = [
    [.40, "#ff4000", True],
    [.50, "#ff8000", False],
    [.70, "#bfff00", False],
    [.75, "#40ff00", True],
    [.80, "#00ff80", False],
    [.85, "#00bfff", False],
    [.90, "#0040ff", False],
    [.95, "#8000ff", False],
    [1.0, "#ff00ff", False],
]

HRR_FILL_COLOR = [
    [.90, 1.0, "#ff0000"],
    [.80, .90, "#dc4806"],
    [.50, .80, "#608d45"],
    [.40, .50, "#069adc"],
    [.0, .40, "#a7adba"]

]

PLOT_HR_LINES = False
PLOT_WIDTH_LINES = False
PLOT_PEAK_POINTS = False


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

    def _load_hr_interval_value(self, _d):
        try:
            self._hr_interval = int(re.search(r'Interval=([0-9]+)', _d).group(1))
        except Exception as e:
            print("Replace this...", e)
            self._invalid = True
            self._invalids.append("Could not load hr_interval")

    def _strip_data(self):
        strip_at = len(self._hr_data)

        try:
            df, df_u = self.get_peaks()

            last_peak = df['peaks'].values[-1]

            for i in range(last_peak, len(self._hr_data)):
                y = self._hr_data[i]
                if self.get_hrr(y) <= (END_HRR_STRIP * 100):
                    strip_at = i
                    break
        except Exception as e:
            print(e)

        self._hr_data = self._hr_data[:strip_at]
        self._load_hr_interval()

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
        self._strip_data()

    def get_max_hr_in_data(self):
        return np.amax(self._hr_data)

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
        try:
            return savgol_filter(self._hr_data, 13, 2)
        except:
            return self._hr_data

    def get_peaks(self):

        def find_bottom(data, peak, reverse=False):
            LOWEST = peak
            iter = range(peak, 0, -1) if reverse else range(peak, len(data))
            for i in iter:
                if data[i] > data[LOWEST]:
                    return LOWEST
                LOWEST = i
            return LOWEST

        def detect_end(data, p_s, p_e, p):
            p_s = p_s.astype(dtype=np.int)
            p_e = p_e.astype(dtype=np.int)
            d_s = np.take(data, p_s)
            d_e = np.take(data, p_e)

            stacked = np.stack((d_s, d_e), axis=0)

            p_amax = np.argmax(stacked, axis=0)
            """
            If value in p_amax is 0, means that the peak_start value is Above the End point
            In this case, we would like to measure width from the equivalent hight on the opposite side of the peak
            """
            new_p_e = np.zeros(shape=p_e.shape, dtype=np.int)
            for i, v in enumerate(p_amax):
                if v == 0:
                    """Find point on other side of the peak."""
                    data_selection = data[p[i]:p_e[i]]

                    lowest_index_local = np.argmin(np.abs(np.subtract(d_s[i], data_selection)))
                    lowest_index_global = p[i] + lowest_index_local
                    new_p_e[i] = lowest_index_global
                else:
                    new_p_e[i] = p_e[i]

            return new_p_e

        def intersect_end(data, peak_start, peak_index):
            for idx in range(peak_index, len(data)):
                v_left = data[peak_start]
                v_right = data[idx]
                if v_right - v_left <= 0:
                    return idx

        if METHOD == 1:
            y = self.smooth_filter()
            low_threshold = self.get_hrr_hr(HRR_PARAM_LOW)
            absolute_lowest = .50
            al_hr = self.get_hrr_hr(absolute_lowest)
            y[y < al_hr] = al_hr

            peaks, _ = find_peaks(y,
                                  height=self.get_hrr_hr(absolute_lowest),
                                  prominence=[12, 50],
                                  width=[
                                      float(PEAK_WIDTH_LOW / self._hr_interval),
                                      float(PEAK_WIDTH_HIGH / self._hr_interval)
                                  ])

            peaks_width, peaks_height, peaks_start, peaks_end = peak_widths(y, peaks, rel_height=0.9)
            peaks_start = peaks_start.astype(np.int)
            peaks_end = peaks_end.astype(np.int)
            peaks_value = np.array([y[p] for p in peaks])
            peaks_valley = np.array([y[p] for p in np.maximum(peaks_start, peaks_end)])

            """ Fix for LARGE LAST PEAK case"""
            if len(peaks_width) > 2:
                if peaks_width[-1] * self._hr_interval > PEAK_WIDTH_HIGH:
                    peaks_start[-1] = peaks_end[-2]

                    old_end = peaks_end[-1]
                    potential_end = intersect_end(y, peaks_start[-1], peaks[-1])

                    w_old = old_end - peaks_start[-1]

                    try:
                        w_new = potential_end - peaks_start[-1]
                    except:
                        w_new = old_end - peaks_start[-1]

                    if w_new < w_old:
                        peaks_end[-1] = potential_end
                    else:
                        peaks_end[-1] = old_end

                    peaks_width[-1] = np.abs(peaks_start[-1] - peaks_end[-1])

            "Fix 2 For large PEAKS"

            def messed_up_peak_detect():

                tdf = pandas.DataFrame({
                    "peaks_start": peaks_start
                })
                tdf = tdf.diff()

                messed_up_indexes = tdf.index[tdf['peaks_start'] < 0].tolist()

                for mid in messed_up_indexes:
                    peaks_start[mid] = peaks_end[mid - 1]
                    old_end = peaks_end[mid]
                    potential_end = intersect_end(y, peaks_start[mid], peaks[mid])

                    w_old = old_end - peaks_start[mid]

                    try:
                        w_new = potential_end - peaks_start[mid]
                    except:
                        w_new = old_end - peaks_start[mid]

                    if w_new < w_old:
                        peaks_end[mid] = potential_end
                    else:
                        peaks_end[mid] = old_end

                    peaks_width[mid] = np.abs(peaks_start[mid] - peaks_end[mid])

            messed_up_peak_detect()

            df = pandas.DataFrame({
                "peaks": peaks,
                "peaks_times": np.array([self._hr_data_times[x] for x in peaks]),
                "peaks_width": (peaks_width * self._hr_interval),
                "peaks_height": peaks_height,
                "peaks_start": peaks_start,
                "peaks_end": peaks_end,
                "peaks_valley": peaks_valley,
                "peaks_value": peaks_value,
                "interval": self._hr_interval
            })
            df_unfiltered = df.copy()

            """
            # Code that allows some cases where the participant has more than 2 peaks with >= HIGH to have ABSOLUTE LOW for the lower bounds. Ref mail.
            """
            try:
                all_higher = df['peaks_value'] > self.get_hrr_hr(HRR_PARAM_HIGH)
                all_higher = pandas.DataFrame({
                    'col1': all_higher
                })
                all_higher["subgroup"] = (all_higher["col1"] != all_higher["col1"].shift(1)).cumsum()
                all_higher["count"] = all_higher.groupby(["subgroup"]).agg('count')
                all_higher = all_higher.dropna()
                meet_criteria = False
                if all_higher.size != 0:

                    last_type = all_higher["col1"].iloc[-1]
                    last_group = all_higher["count"].iloc[-1]

                    if last_type == True and last_group >= 2:
                        meet_criteria = True

                if meet_criteria:
                    required_lower_bound = absolute_lowest
                else:
                    required_lower_bound = HRR_PARAM_LOW
            except:
                required_lower_bound = HRR_PARAM_LOW

            """Width constraint."""
            df = df[(df['peaks_width'] >= PEAK_WIDTH_LOW) & (df['peaks_width'] <= PEAK_WIDTH_HIGH)]

            """Height Threshold"""
            df = df[(df['peaks_value']) >= self.get_hrr_hr(required_lower_bound)]

            """Height constraint."""
            df = df[(df['peaks_height'] >= PEAK_HEIGHT)]

            return df, df_unfiltered

        elif METHOD == 2:

            x = self.smooth_filter()
            import peakutils

            peaks = np.array(peakutils.indexes(x))
            peaks_start = np.array([find_bottom(x, peak, reverse=True) for peak in peaks])
            peaks_end = np.array([find_bottom(x, peak) for peak in peaks])
            peaks_end = detect_end(x, peaks_start, peaks_end, peaks)  # Refine end for some cases
            peaks_value = np.array([x[p] for p in peaks])

            peaks_width = np.abs(peaks_start - peaks_end)
            peaks_valley = np.array([x[p] for p in np.maximum(peaks_start, peaks_end)])

            peaks_height = np.abs(peaks_value - peaks_valley)

            df = pandas.DataFrame({
                "peaks": peaks,
                "peaks_times": np.array([self._hr_data_times[x] for x in peaks]),
                "peaks_width": (peaks_width * self._hr_interval),
                "peaks_height": peaks_height,
                "peaks_start": peaks_start,
                "peaks_end": peaks_end,
                "peaks_valley": peaks_valley,
                "peaks_value": peaks_value,
                "interval": self._hr_interval
            })
            df_unfiltered = df.copy()

            """Width constraint."""
            df = df[(df['peaks_width'] >= PEAK_WIDTH_LOW) & (df['peaks_width'] <= PEAK_WIDTH_HIGH)]
            """Height Threshold"""
            df = df[(df['peaks_value']) >= self.get_hrr_hr(HRR_PARAM_LOW)]

            """Height constraint."""
            df = df[(df['peaks_height'] >= PEAK_HEIGHT)]

            """if "16060201" in self._path:
    
                print(df_unfiltered)
                print(df)
    
                #y = df["peaks"].values
                #x = np.arange(0, len(y))
    
                plt.cla()
                plt.clf()
                plt.plot(x)
                plt.plot(df["peaks"].values, [x[_] for _ in df["peaks"].values], "s", color="#ffa500")
                for p in range(len(df["peaks"].values)):
                    plt.hlines(max(x[df['peaks_start'].values[p]], x[df['peaks_end'].values[p]]), xmin=df['peaks_start'].values[p], xmax=df['peaks_end'].values[p], color="C1")
                plt.savefig("test.png")
    
                print(len(df["peaks"].values))
    
            _ = {}"""

            return df, df_unfiltered

    def get_values_over_hrr(self, hrr=0):
        bpm = self.get_hrr_hr(hrr)
        data = np.copy(self._hr_data)
        data[self._hr_data < bpm] = 0
        return data

    def plot(self, ax=None, save=False, save_location=None, save_name=None, no_labels=False):
        if save is True and (save_location is None or save_name is None):
            raise RuntimeError("You must set a save location AND save_name")

        fig = None
        if ax is None:
            fig, ax = plt.subplots(1, 1)

        hr_data = self._hr_data  # self.smooth_filter()
        hr_data_times = self._hr_data_times
        hr_interval = self._hr_interval
        df, df_u = self.get_peaks()

        hr_data_times *= hr_interval
        hr_data_times = np.divide(hr_data_times, 60)

        """Plot background."""
        for spec in HRR_FILL_COLOR:
            _low = self.get_hrr_hr(spec[0])
            _high = self.get_hrr_hr(spec[1])
            _color = spec[2]
            plt.fill_between(hr_data_times, _low, _high, color=_color, alpha=0.5)

        """Plot lines."""
        sns.lineplot(hr_data_times, hr_data, color="red", linestyle="-", linewidth=1.0)

        """Fill under line"""
        ax.fill_between(hr_data_times, hr_data, self.get_hrr_hr(.0), interpolate=True, color='#E0FFFF', alpha=0.5)

        """Plot peaks"""
        if PLOT_PEAK_POINTS and no_labels is False:
            for p in df["peaks"]:
                is_low = self.get_hrr(hr_data[p]) < (HRR_PARAM_HIGH * 100)

                if is_low:
                    ax.plot(hr_data_times[p], hr_data[p], "s", color="#ffa500")
                else:
                    ax.plot(hr_data_times[p], hr_data[p], "s", color="#ff0000")

        """Scale width of hlines properly"""
        """Plot width of peak"""
        if PLOT_WIDTH_LINES and no_labels is False:
            for p in range(len(df["peaks"].values)):
                plt.hlines(
                    max(
                        hr_data[df['peaks_start'].values[p]],
                        hr_data[df['peaks_end'].values[p]]
                    ),
                    xmin=np.multiply(df['peaks_start'].values[p], hr_interval) / 60,
                    xmax=np.multiply(df['peaks_end'].values[p], hr_interval) / 60,
                    color="C1"
                )

        """Plot horizontal marker lines"""
        for _hrrs in HRR_LINES:
            if (PLOT_HR_LINES is False or no_labels is True) and _hrrs[2] is False:
                continue
            plt.axhline(self.get_hrr_hr(_hrrs[0]), color=_hrrs[1], linewidth=1.0,
                        label="%s%s" % (int(_hrrs[0] * 100), "% HRR"), linestyle="--")

            # plt.axhline(self.get_hrr_hr(HRR_PARAM_LOW), color="C2", linewidth=3.0, label="LOW-%s" % HRR_PARAM_LOW, linestyle="--")
            # plt.axhline(self.get_hrr_hr(HRR_PARAM_HIGH), color="C3", linewidth=3.0, label="HIGH-%s" % HRR_PARAM_HIGH, linestyle="--")

        ax.set(xlabel='Duration in Minutes', ylabel='Heart Rate (bps)')
        ax.margins(0.0)

        # start, end = ax.get_ylim()
        # ax.yaxis.set_ticks(np.arange(start, end, 10))

        # plt.show(block=False)

        handles, labels = ax.get_legend_handles_labels()
        handles.insert(0, mpatches.Patch(color='none', label="HRmax: %s" % self._max_hr))

        plt.legend(handles=handles, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

        if PLOT_PEAK_POINTS and no_labels is False:
            plt.title(os.path.basename(self._path))
        plt.tight_layout()
        plt.savefig(os.path.join(save_location, save_name), bbox_inches="tight", pad_inches=0, dpi=DPI)
        if SAVE_PDF:
            plt.savefig(os.path.join(save_location, save_name) + ".pdf", bbox_inches="tight", pad_inches=0, dpi=DPI)
        plt.close(fig)

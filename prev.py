import os
import glob
import zipfile
import ntpath
import shutil
import re
import logging
import datetime
from collections import OrderedDict

import pandas as pd
import math
pd.set_option('display.expand_frame_repr', False)
import numpy as np
import matplotlib.pyplot as plt
import six

mpl_logger = logging.getLogger("matplotlib")
mpl_logger.setLevel(logging.WARNING)

logging.basicConfig(level=logging.DEBUG)
_LOGGER = logging.getLogger("cancer_parser")


dir_path = os.path.dirname(os.path.realpath(__file__))
ROLLING_AVERAGE = False



BREATHING_COLUMN = ["RF", "BF", "RR"]
EXCEL_EXTENSIONS = (".xlsx", ".xls")
EXCEL_DESTINATION = "excel"
EXCEL_TEMP_DESTINATION = "excel_tmp"
CSV_DESTINATION = "csv"
ZIP_DESTINATION = "zips"
FIGURE_DESTINATION = "figures"
CUT_OFF_THRESHOLD = .20
SAMPLE_DURATION = 60.0
PLATEU_CRITERIA = {
    "150 ml/min": lambda x: 1 if x <= .15 else 0,
    "80 ml/min": lambda x: 1 if x <= .08 else 0,
    "50 ml/min": lambda x: 1 if x <= .05 else 0
}

EARLY_PLATEU = 0.7
PLATEU_SAMPLE_SIZE_BEFORE = 1
PLATEU_SAMPLE_SIZE_AFTER = 1
INVALID_COLUMN_THRESHOLD = .8  # This threshold is used to remove useless stuff that is placed in columns before the table
INVALID_ROW_THRESHOLD = .40  # This threshold is used to remove string / header  in first clean iteration
INVALID_ROW_TYPE_THRESHOLD = .50  # Used to determine if a row contains invalid data. such as not numbers in the VO2 column


class Runtime:

    def __init__(self):
        self.total_missing_t3 = 0

    def ensure_paths(self):
        os.makedirs(os.path.join(dir_path, CSV_DESTINATION), exist_ok=True)
        os.makedirs(os.path.join(dir_path, ZIP_DESTINATION), exist_ok=True)
        os.makedirs(os.path.join(dir_path, FIGURE_DESTINATION), exist_ok=True)
        os.makedirs(os.path.join(dir_path, EXCEL_DESTINATION), exist_ok=True)

    @staticmethod
    def render_mpl_table(data, col_width=3.0, row_height=0.625, font_size=14,
                         header_color='#40466e', row_colors=['#f1f1f2', 'w'], edge_color='w',
                         bbox=[0, 0, 1, 1], header_columns=0,
                         ax=None, **kwargs):
        if ax is None:
            size = (np.array(data.shape[::-1]) + np.array([0, 1])) * np.array([col_width, row_height])
            fig, ax = plt.subplots(figsize=size)
            ax.axis('off')

        mpl_table = ax.table(cellText=data.values, bbox=bbox, colLabels=data.columns, **kwargs)

        mpl_table.auto_set_font_size(False)
        mpl_table.set_fontsize(font_size)

        for k, cell in six.iteritems(mpl_table._cells):
            cell.set_edgecolor(edge_color)
            if k[0] == 0 or k[1] < header_columns:
                cell.set_text_props(weight='bold', color='w')
                cell.set_facecolor(header_color)
            else:
                cell.set_facecolor(row_colors[k[0]%len(row_colors) ])
        return ax

    def unzipper(self, file_path=os.path.join(dir_path, ZIP_DESTINATION)):

        for zip_file_path in glob.glob(os.path.join(file_path, "*.zip")):

            zip_ref = zipfile.ZipFile(zip_file_path, 'r')

            excel_files = [x for x in zip_ref.filelist if x.filename.endswith(EXCEL_EXTENSIONS)]

            for excel_file in excel_files:
                excel_file_dest = os.path.join(dir_path, EXCEL_DESTINATION, ntpath.basename(excel_file.filename))
                with zipfile.ZipFile(zip_file_path) as z:
                    with open(excel_file_dest, 'wb') as f:
                        f.write(z.read(excel_file.filename))
            zip_ref.close()

    def generate_patients(self):
        errors = {
            'patient_id': [],
            't_id': []
        }
        patient_map = {}

        for f in sorted([x for x in os.listdir(os.path.join(dir_path, EXCEL_DESTINATION)) if x.endswith(EXCEL_EXTENSIONS)]):
            patient_id = re.findall('^[a-zA-Z0-9]{5}', f)

            if not patient_id:
                errors["patient_id"].append(f)
                continue

            patient_id = patient_id[0]

            if patient_id not in patient_map:
                patient_map[patient_id] = Patient(patient_id, self)

            t_id = re.findall('([[t,T][0-4]{1}|000[0-4]{1})', f)

            if not t_id:
                patient_map[patient_id].errors_t_id.append(f)
                errors["t_id"].append(f)
                continue

            patient_map[patient_id].add_t_id(t_id, f)

        _LOGGER.debug("Patient ID errors: %s" % errors["patient_id"])
        _LOGGER.debug("T_ID errors: %s" % errors["t_id"])
        _LOGGER.debug("Total patient file errors: %s" % (len(errors["t_id"]) + len(errors["patient_id"])))
        return patient_map


class Patient:

    def __init__(self, id, runtime):
        self.id = id
        self.runtime = runtime

        self.data = {
            "T0": [],
            "T1": [],
            "T2": [],
            "T3": [],
            "T4": []
        }

        self.results = {
            "T0": [],
            "T3": [],
            "T4": []
        }

        self.errors_t_id = []

    def plot(self):
        merged_results = self.results["T0"] + self.results["T3"] + self.results["T4"]
        fig, axes = plt.subplots(nrows=len(merged_results)*2, ncols=1)
        fig.set_figheight(20)
        fig.set_figwidth(10)

        for i, data in enumerate(merged_results):



            if isinstance(axes, np.ndarray):
                ax = data["data"]["vo2"].plot(ax=axes[i])
            else:
                ax = data["data"]["vo2"].plot()
            ax.set_title(data["title"])
            ax.set_xlabel('Minutes')
            ax.set_ylabel('VO2')
            ax.set_xticks(np.arange(0, len(data["data"])+1, 1.0))
            ax.legend()
            ax.grid()



            colors = ["r", "g", "b"]

            ymin, ymax = ax.get_ylim()
            for c, plateu_type in enumerate(reversed(list(PLATEU_CRITERIA.keys()))):

                if (data["data"][plateu_type] == 1).any():
                    # If plateu is reached for this type
                    idx = data["data"][plateu_type].idxmax()
                    color = colors[c]

                    ax.vlines(x=[idx], ymin=ymin, ymax=ymax, color=color, label=plateu_type)


        for i, data in enumerate(merged_results):
            Runtime.render_mpl_table(data["data"], ax=axes[i+len(merged_results)])

        plt.tight_layout()
        fig.savefig(os.path.join(dir_path, "figures", self.id + ".png"))#, dpi=800)#, dpi=1000)
        plt.close()

    def add_t_id(self, t_id, file):

        t_id_list = list(reversed(sorted(t_id)))

        # Prioritize capital id
        if len(t_id_list) == 2:

            if t_id_list[0][0].islower() and t_id_list[1][0].isupper():
                t_id = t_id_list[1]
            else:
                t_id = t_id_list[0]
        elif len(t_id_list) == 1:
            t_id = t_id_list[0].upper()

            if t_id == "0001":
                t_id = "T1"
            elif t_id == "0002":
                t_id = "T3"
            #else:
            #    raise ValueError("Should be 0001 or 0002. got: %s" % t_id)

        self.data[t_id].append(file)

    def _find_correct_file(self, lst):
        found_files = 0
        data_frames = []
        for it in lst:
            file_path = os.path.join(dir_path, EXCEL_DESTINATION, it)
            df = pd.read_excel(file_path, header=None)

            data = df.values  # Convert to numpy

            # Clean useless rows
            rows_to_remove = []
            for row_i, row in enumerate(data):
                num_nan = pd.isnull(row)  # Find number of empty cells for row
                nan_percent = list(num_nan).count(True) / len(num_nan)

                if nan_percent >= INVALID_ROW_THRESHOLD:
                    rows_to_remove.append(row_i)

            """# Clean useless columns
            columns_to_remove = []
            for col_i, column in enumerate(data.T):
                num_nan = pd.isnull(column)  # Find number of empty cells for row
                nan_percent = list(num_nan).count(True) / len(num_nan)

                if nan_percent >= INVALID_COLUMN_THRESHOLD:
                    columns_to_remove.append(col_i)


            data = np.delete(data, columns_to_remove, 1)"""
            data = np.delete(data, rows_to_remove, 0)

            # Clean / Strip spaces from string based cells
            data = np.array([[item.strip() if isinstance(item, str) else item for item in s ] for s in data])

            # Find v02 locations
            vo2_locations = []
            for n in ["V'O2", "VO2"]:
                for x, y in zip(*np.where(data == n)):
                    vo2_locations.append((x, y))

            #, ensure that all vo2locations are at the same y position?

            if not vo2_locations and "puls" not in file_path.lower():
                _LOGGER.debug("The file is invalid. Please check:  %s", file_path)
                continue

            if "puls" in file_path.lower():
                continue

            # Skew check
            vo2_locations_y = [el[1] for el in vo2_locations]
            if vo2_locations_y.count(vo2_locations_y[0]) != len(vo2_locations_y):
                _LOGGER.debug("The V02 field is skewed which makes this kinda hard...")
                continue

            # Remove string values
            rows_to_remove = []
            for row_i, row in enumerate(data):
                cast_count = 0  # Number of successful casts
                for cell in row:

                    # Test if number
                    try:
                        float(cell)
                        cast_count += 1
                        continue
                    except:
                        pass

                    if isinstance(cell, float) or isinstance(cell, int):
                        cast_count += 1
                        continue

                    if isinstance(cell, datetime.time) or isinstance(cell, datetime.datetime) or ":" in cell:
                        # Probably time
                        cast_count += 1
                        continue

                cast_percent = cast_count / len(row)

                if cast_percent < INVALID_ROW_TYPE_THRESHOLD:
                    if row_i == 0:  # TODO skip header?
                        continue
                    rows_to_remove.append(row_i)

            data = np.delete(data, rows_to_remove, 0)

            df = pd.DataFrame(data[1:], columns=data[0])

            found_files += 1

            data_frames.append((file_path, df))
        return data_frames

    @staticmethod
    def time_convert(x):
        if isinstance(x, datetime.time):

            return ((x.hour * 60) * 60) + (x.minute * 60) + x.second
        elif isinstance(x, str):
            times = x.split(":")
            if len(times) == 2:
                return int((int(times[0]) * 60) + int(times[1]))
            elif len(times) == 3:
                return int((int(times[0]) * 60 * 60) + (int(times[1]) * 60) + int(times[2]))
            else:
                return 0  # TODO handling

        else:
            return 0 # TODO handling

    def calculate(self, test_type, filename, data_frame):
        results = OrderedDict()
        results.update({
            "id": self.id,
            "stage": test_type
        })
        results.update({
            k: "N/A" for k in PLATEU_CRITERIA.keys()
        })



        #####################################
        # Remove duplicate columns (ie two Time columns)
        #####################################
        data_frame = data_frame.loc[:,~data_frame.columns.duplicated()]

        #####################################
        # Identify column names
        #####################################
        time_column_name = [t for t in data_frame.columns if 'time' in str(t).lower()][0]
        vo2_column_name = [s for s in data_frame.columns if any(xs == s for xs in ["V'O2", "VO2"])][0]

        #####################################
        # Retrieve vo2 column records
        #####################################
        vo2_df = data_frame[vo2_column_name].replace(r'-', 0, regex=True).astype(float).rename("vo2")

        # Convert from vo2/ml to vo2/l
        vo2_df = vo2_df.divide(1000) if vo2_df.mean() > 50 else vo2_df

        #####################################
        # Retrieve time column records
        #####################################
        duration_df = data_frame[time_column_name].apply(self.time_convert).rename("duration")

        #####################################
        # Compute cut-off
        #####################################

        vo2_pct_change_df = vo2_df.pct_change().abs().fillna(0).rename("cut-off")

        cut_off_pd = pd.concat([vo2_df, duration_df, vo2_pct_change_df], axis=1)
        cut_off_pd = cut_off_pd[cut_off_pd["cut-off"] <= CUT_OFF_THRESHOLD]

        # Reallign vo2_df and duration_df with filtered objects
        vo2_df = cut_off_pd["vo2"]
        duration_df = cut_off_pd["duration"]

        #####################################
        # Compute difference between timesteps
        #####################################
        diff_df = duration_df.diff(periods=1).fillna(0).rename("time_diff")
        dominating_value = int(diff_df.mode().max())

        #####################################
        # Generate a group id, that groups each record into a minute
        #####################################
        time_group_df = duration_df.divide(60).apply(np.ceil).apply(lambda x: 1 if x == 0 else x).astype(int).rename("group_id")

        #####################################
        # Calculate the vo2 mean over the group_id
        #####################################
        vo2_mean_min_df = pd.concat([vo2_df, time_group_df], axis=1).groupby(["group_id"]).mean()

        #####################################
        # Calculate the minutes sum for each group
        #####################################
        duration_sum_min_df = pd.concat([diff_df, time_group_df], axis=1).groupby(["group_id"]).sum()

        #####################################
        # Find the maximum vo2 position
        #####################################
        if ROLLING_AVERAGE:
            #####################################
            # Calculate the vo2 mean and time in the rolling window
            #####################################
            rolling_means = []
            rolling_means_start = []
            rolling_means_end = []
            total_time = 0
            for i in range(len(vo2_df.values)):
                cum_time = 0
                mean_vals = []
                for j in range(len(vo2_df.values)):
                    if i + j >= len(diff_df.values):
                        break

                    t = diff_df.values[i + j]
                    v = vo2_df.values[i + j]

                    mean_vals.append(v)
                    cum_time += t

                    if cum_time >= SAMPLE_DURATION:
                        break

                rolling_means.append(np.array(mean_vals).mean())
                rolling_means_start.append(total_time)
                rolling_means_end.append(total_time + cum_time)
                total_time += cum_time

            rolling_means = pd.DataFrame(data={
                'vo2': rolling_means,
                'Start': rolling_means_start,
                'End': rolling_means_end,
                'Diff': np.subtract(rolling_means_end, rolling_means_start)
            })

            # Strip away times that does not conform the time-scale  # TODO might be wrong?

            rolling_means = rolling_means[rolling_means["Diff"] == SAMPLE_DURATION]

            vo2_check_df = rolling_means["vo2"]
            duration_check_df = rolling_means["Diff"]
        else:
            vo2_check_df = vo2_mean_min_df["vo2"]
            duration_check_df = duration_sum_min_df

        vo2_check_df = vo2_check_df.reset_index(drop=True)

        # Ensure that there is any records left to check.
        if len(vo2_check_df) == 0:
            _LOGGER.info("NO Records to use")
            return results

        max_value_idx = int(vo2_check_df.idxmax())


        early_detect = "Yes" if max_value_idx / len(vo2_check_df) < EARLY_PLATEU else "No"
        check_scope = vo2_check_df.iloc[max_value_idx - 1 : max_value_idx + 2]
        check_scope = pd.concat([check_scope, check_scope.diff().abs().rename("diff")], axis=1)

        for plateu_type, criterion in PLATEU_CRITERIA.items():
            check_scope = pd.concat([
                check_scope,
                check_scope["diff"].apply(criterion).rename(plateu_type)
            ], axis=1)

        # Check if each plateu column has a 1 value
        for plateu_type, criterion in PLATEU_CRITERIA.items():
            results[plateu_type] = 1 if (check_scope[plateu_type] == 1).any() else 0

        results["id"] = self.id
        results["stage"] = test_type
        results["early"] = early_detect

        results["title"] = self.id
        results["data"] = pd.concat([check_scope.drop("vo2", axis=1).drop("diff", axis=1), vo2_mean_min_df], axis=1)
        return results

    def get_breathing(self, test_type):
        if test_type == "T0":
            dfs = self.t1_dfs
        if test_type == "T3":
            dfs = self.t3_dfs

        accum = []
        for test in dfs:
            test_filename = test[0]
            test_df = test[1]

            try:
                breathing_column_name = [t for t in test_df.columns if str(t).upper() in BREATHING_COLUMN][0]
            except IndexError:
                # Could not find breathing
                return [self.id, "N/A"]

            breathing_data = test_df[breathing_column_name]

            data = [self.id] + breathing_data.values.tolist()
            accum.extend(data)

        return accum


    def find(self):
        self.t1_dfs = self._find_correct_file(self.data["T0"])
        self.t3_dfs = self._find_correct_file(self.data["T3"])
        self.t4_dfs = self._find_correct_file(self.data["T4"])

    def inspect(self):
        t1_dfs = self._find_correct_file(self.data["T0"])
        t3_dfs = self._find_correct_file(self.data["T3"])
        t4_dfs = self._find_correct_file(self.data["T4"])
        self.results["T0"] = [self.calculate("T0", *_) for _ in t1_dfs]
        self.results["T3"] = [self.calculate("T3", *_) for _ in t3_dfs]
        self.results["T4"] = [self.calculate("T4", *_) for _ in t4_dfs]

        #[self.get_breathing("T0", *_) for _ in t1_dfs]
        #[self.get_breathing("T3", *_) for _ in t3_dfs]
        #[self.get_breathing("T4", *_) for _ in t4_dfs]

        """self.results["T0"] = [self.calculate("T0", *_) for _ in self.t1_dfs]
        self.results["T3"] = [self.calculate("T3", *_) for _ in self.t3_dfs]
        self.results["T4"] = [self.calculate("T4", *_) for _ in self.t4_dfs]"""



if __name__ == "__main__":

    # Create runtime intance that keeps track of statistics. etc
    runtime = Runtime()

    runtime.ensure_paths()
    runtime.unzipper()
    patient_map_generator = runtime.generate_patients()

    _LOGGER.debug("Found %s patients.",  len(patient_map_generator))


    accumulated_results = []
    for (patient_id, patient) in patient_map_generator.items():
        patient.inspect()
        patient.plot()

        for v in patient.results.values():
            accumulated_results.extend(v)

    writer = pd.ExcelWriter('output.xlsx')
    df = pd.DataFrame(accumulated_results)
    df.to_excel(writer, "Results")
    writer.save()

    """
    accumulated_results = []
    for (patient_id, patient) in patient_map_generator.items():
        patient.find()

        data = patient.get_breathing("T0")
        if data:
            accumulated_results.append(data)

    writer = pd.ExcelWriter('breathing.xlsx')
    df = pd.DataFrame(accumulated_results)
    print(df)

    df.to_excel(writer, "Results")
    writer.save()
    """
    print("TOTAL MISSING T3: %s " % runtime.total_missing_t3)

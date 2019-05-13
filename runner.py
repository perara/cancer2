import collections
import glob
import io
import multiprocessing
import random
import shutil
import zipfile
import re
import logging
from random import shuffle
import xlsxwriter
import pandas as pd
import os
import numpy as np
from util import HRM
import util
from datetime import datetime
# 2 UKER - 1 MAI
dir_path = os.path.dirname(os.path.realpath(__file__))

pd.set_option('display.expand_frame_repr', False)

PATH_MUST_CONTAIN = ["Hög", "HÖG"]
HR_FILE_EXTENSIONS = ('.hrm')
LABELING_MODE = True
LABELING_N = 100

ZIP_DESTINATION = os.path.join(dir_path, "zips")
OUTPUT_LOCATION = os.path.join(dir_path, "output")
HRM_DESTINATION = os.path.join(OUTPUT_LOCATION, "HRM")
EXCEL_DEST = os.path.join(OUTPUT_LOCATION, "EXCEL")
PLOT_DESTINATION = os.path.join(OUTPUT_LOCATION, "PLOTS")
PLOT_NO_LABELS_DESTINATION = os.path.join(OUTPUT_LOCATION, "PLOTS_NO_LABELS")
LABELING_DESTINATION = os.path.join(OUTPUT_LOCATION, "LABELING")
PLOT_VALID_DESTINATION = os.path.join(PLOT_DESTINATION, "valid")
PLOT_INVALID_DESTINATION = os.path.join(PLOT_DESTINATION, "invalid")

CONSTRAINT_DURATION = 15 * 60
CONSTRAINT_N_PEAKS = 5

os.makedirs(OUTPUT_LOCATION, exist_ok=True)
os.makedirs(HRM_DESTINATION, exist_ok=True)
os.makedirs(EXCEL_DEST, exist_ok=True)
os.makedirs(PLOT_DESTINATION, exist_ok=True)
os.makedirs(PLOT_NO_LABELS_DESTINATION, exist_ok=True)
os.makedirs(LABELING_DESTINATION, exist_ok=True)

shutil.rmtree(PLOT_INVALID_DESTINATION, ignore_errors=True)
shutil.rmtree(PLOT_VALID_DESTINATION, ignore_errors=True)
shutil.rmtree(LABELING_DESTINATION, ignore_errors=True)

os.makedirs(PLOT_INVALID_DESTINATION, exist_ok=True)
os.makedirs(PLOT_VALID_DESTINATION, exist_ok=True)
os.makedirs(LABELING_DESTINATION, exist_ok=True)

mpl_logger = logging.getLogger("matplotlib")
mpl_logger.setLevel(logging.WARNING)

logging.basicConfig(level=logging.DEBUG)
_LOGGER = logging.getLogger("cancer_parser")


def get_patient_id(haystack):
    rex_1 = re.search(r'[0-9]{3}[a-zA-Z]{2}', haystack)
    rex_2 = re.search(r'[0-9]{3} [a-zA-Z]{2}', haystack)
    rex_3 = re.search(r'[0-9]{2}[a-zA-Z]{3}', haystack)

    # Special case: Anita Modigh
    rex_4 = re.search(r'Anita Modigh', haystack)

    if rex_1:
        return rex_1.group(0)
    elif rex_2:
        return rex_2.group(0).replace(" ", "")
    elif rex_3:
        return rex_3.group(0)
    elif rex_4:
        return "000AA"  # TODO note this
    else:
        raise NotImplementedError("A case were not implemented: %s" % haystack)


def get_year(haystack):
    rex_1 = re.search(r'20[0-9]{2}', haystack)

    if rex_1:
        return rex_1.group(0)
    else:
        return 0000


def get_hrm_data(zip_file_path, hrm_filepath):
    with zipfile.ZipFile(zip_file_path) as z:
        return z.read(hrm_filepath)


def save_data(hrm_id, hrm_year, hrm_user, hrm_data):
    hrm_filename = "%s-%s-%s.hrm" % (hrm_user, hrm_year, hrm_id)
    hrm_full_path = os.path.join(HRM_DESTINATION, hrm_filename)

    with open(hrm_full_path, 'wb') as f:
        f.write(hrm_data)


def unzipper(file_path=ZIP_DESTINATION, ignore=True):
    l = list(glob.glob(os.path.join(file_path, "*.zip")))
    shuffle(l)
    for zip_file_path in l:

        # Open ZIP
        zip_ref = zipfile.ZipFile(zip_file_path, 'r')

        # Find all hrm files (absolute path)
        hrm_files = [x for x in zip_ref.filelist if
                     x.filename.endswith(HR_FILE_EXTENSIONS) and any(s in x.filename for s in PATH_MUST_CONTAIN)]
        hrm_dest_n_files = len(list(os.listdir(HRM_DESTINATION)))

        if len(hrm_files) == hrm_dest_n_files or ignore:
            """Do not reparse everything..."""
            return True

        for hrm_object in hrm_files:
            hrm_filepath = hrm_object.filename
            file_name_chunks = hrm_filepath.split("/")

            hrm_id = file_name_chunks[-1].split(".")[0].split(" ")[0]
            hrm_year = get_year(hrm_filepath)
            hrm_user = get_patient_id(hrm_filepath)
            hrm_data = get_hrm_data(zip_file_path, hrm_filepath)

            save_data(hrm_id, hrm_year, hrm_user, hrm_data)

        zip_ref.close()


class Stats:

    def __init__(self):
        self.average_peak_widths = []

    def update_avg_peak_widths(self, a, b, peak_widths):
        self.average_peak_widths.extend(peak_widths[0])


def worker(hrm_filepath):
    # print(hrm_filepath)
    file_name = os.path.basename(hrm_filepath)
    file_name = os.path.splitext(file_name)[0]
    splt = file_name.split("-")
    hrm_user = splt[0]
    hrm_year = splt[1]
    hrm_id = splt[2]

    hrm = HRM(hrm_filepath)
    if hrm.validate():
        hrm.plot(save=True, save_location=PLOT_DESTINATION, save_name=file_name + ".png")

        if LABELING_MODE:
            hrm.plot(save=True, save_location=PLOT_NO_LABELS_DESTINATION, save_name=file_name + ".png", no_labels=True)

        df, df_u = hrm.get_peaks()
        hrm_peak_times = df["peaks_times"].values
        hrm_peak_types = [
            "LOW" if hrm.get_hrr(x) < (util.HRR_PARAM_HIGH * 100) else "HIGH" for x in df["peaks"]
        ]

        if len(df["peaks_end"]) > 0:

            first_end = df["peaks_end"].values[0]
            last_end = df["peaks_end"].values[-1]
            interval_data = hrm._hr_data[first_end:last_end]
        else:
            interval_data = np.array([])

        return dict(
            user=hrm_user,
            year=hrm_year,
            id=hrm_id,
            peak_widths=df["peaks_width"].values,
            peaks_values=df["peaks_value"].values,
            peaks_lowpoints=df["peaks_valley"].values,
            peaks_interval_data=interval_data,
            peaks_session=np.trim_zeros(hrm.get_values_over_hrr(.40)),
            interval_scale=hrm._hr_interval,
            peak_types=hrm_peak_types,
            max_hr=hrm.get_max_hr_in_data(),
            duration=hrm.get_duration(),
            n_peaks=len(df["peaks"].values),
            peak_times=hrm_peak_times,
            file_name=file_name,
            date=datetime.strptime(hrm_id[0:6], "%y%m%d"),
            df=df,
            df_u=df_u

        )

    return None



def insert_into_workbook(workbook, worksheet, x, res, valid, images, reasons=[]):
    bold = workbook.add_format({'bold': True})
    date_format = workbook.add_format({'num_format': 'm/dd/yyyy'})
    def print_pandas(x_start, y_start, df, sheet):
        cols = df.columns.values

        for i, col in enumerate(cols):
            sheet.write(y_start, x_start + i, col, bold)

        for index, row in df.iterrows():
            values = row.values
            for i, val in enumerate(values):
                sheet.write(y_start + 1 + index, x_start + i, val)

    if images:
        sheet_name = str(res["id"]) + str(random.getrandbits(15))
        images_sheet = workbook.add_worksheet(name=sheet_name)
        worksheet.write(x, 0, "=HYPERLINK(\"#%s!B1\", \"Chart\")" % sheet_name)
        plot_dest = PLOT_VALID_DESTINATION if valid else PLOT_INVALID_DESTINATION
        images_sheet.write(0, 1, "=HYPERLINK(\"#Main!K%s\", \"GO BACK\")" % (x + 1))
        images_sheet.write(1, 1, "Failure Reason", bold)
        images_sheet.write(1, 2, ','.join(reasons))

        images_sheet.insert_image("B5", os.path.join(plot_dest, res["file_name"] + ".png"))

        print_pandas(10, 2, res["df"], images_sheet)
        print_pandas(10, 2 + 1 + len(res["df"]), res["df_u"], images_sheet)

    worksheet.write(x, 1, res["user"])
    worksheet.write(x, 2, str(res["id"]))
    worksheet.write(x, 3, str(res["year"]))
    worksheet.write(x, 4, res["file_name"])
    worksheet.write(x, 5, (str(res["duration"] / 60)))
    worksheet.write(x, 6, (str(res["n_peaks"])))
    worksheet.write(x, 7, str(res["peaks_values"].mean()).replace("nan", "-"))
    worksheet.write(x, 8, str(res["peaks_lowpoints"].mean()).replace("nan", "-"))
    worksheet.write(x, 9, str(res["peaks_session"].mean()).replace("nan", "-"))
    worksheet.write(x, 10, str(res["peaks_interval_data"].mean()).replace("nan", "-"))
    worksheet.write(x, 11, str(res["peak_widths"].mean()).replace("nan", "-"))
    worksheet.write(x, 12, str(res["max_hr"]).replace("nan", "-"))
    worksheet.write(x, 13, int(valid))
    worksheet.write_datetime(x, 14, res["date"], date_format)
    worksheet.write(x, 15, ', '.join([str(x) for x in np.divide(res["peak_times"], 60)]))


def parse(n_start=0, n_files=-1, images=True, override=None, prefix="output", sort=True):
    stats = Stats()

    """Retrieve files"""
    if override is not None:
        files = override
    else:
        files = list(glob.glob(os.path.join(HRM_DESTINATION, "*.hrm")))
        shuffle(files)
        files = files[n_start:n_files]

    if n_start == 0 and n_files == -1:
        xls_name = os.path.join(EXCEL_DEST, '%s.xlsx') % (prefix)
        xls_session_name = os.path.join(EXCEL_DEST, '%s_mean.xlsx') % (prefix)
        xls_session_a_name = os.path.join(EXCEL_DEST, '%s_mean_week.xlsx') % (prefix)
        xls_session_s_name = os.path.join(EXCEL_DEST, '%s_mean_week_success.xlsx') % (prefix)
    else:
        xls_name = os.path.join(EXCEL_DEST, '%s-%s-%s.xlsx') % (prefix, n_start, n_files)
        xls_session_name = os.path.join(EXCEL_DEST, '%s-%s-%s_mean.xlsx') % (prefix, n_start, n_files)
        xls_session_a_name = os.path.join(EXCEL_DEST, '%s-%s-%s_mean_week.xlsx') % (prefix, n_start, n_files)
        xls_session_s_name = os.path.join(EXCEL_DEST, '%s-%s-%s_mean_week_success.xlsx') % (prefix, n_start, n_files)

    workbook = xlsxwriter.Workbook(xls_name)
    worksheet = workbook.add_worksheet(name="Main")

    workbook_session = xlsxwriter.Workbook(xls_session_name)
    worksheet_session = workbook_session.add_worksheet(name="Main")
    date_format_s = workbook_session.add_format({'num_format': 'm/dd/yyyy'})
    bold_s = workbook_session.add_format({'bold': True})
    worksheet_session.write(0, 0, "User", bold_s)
    worksheet_session.write(0, 1, "Filename", bold_s)
    worksheet_session.write(0, 2, "Date", bold_s)
    worksheet_session.write(0, 3, "Peak Count", bold_s)
    worksheet_session.write(0, 4, "Peak Heights(Mean)", bold_s)
    worksheet_session.write(0, 5, "Lowpoint (Mean)", bold_s)
    worksheet_session.write(0, 6, "Session > 40% (Mean)", bold_s)
    worksheet_session.write(0, 7, "Interval (Mean)", bold_s)
    worksheet_session.write(0, 8, "Peak Widths (Mean)", bold_s)
    worksheet_session.write(0, 9, "Accepted", bold_s)
    worksheet_session.write(0, 10, "Year", bold_s)

    workbook_session_a = xlsxwriter.Workbook(xls_session_a_name)
    worksheet_session_a = workbook_session_a.add_worksheet(name="Main")
    bold_a = workbook_session_a.add_format({'bold': True})
    worksheet_session_a.write(0, 0, "User", bold_a)
    worksheet_session_a.write(0, 1, "Session Count", bold_a)
    worksheet_session_a.write(0, 2, "Duration (Sum)", bold_a)
    worksheet_session_a.write(0, 3, "Duration (Mean)", bold_a)
    worksheet_session_a.write(0, 4, "Peaks (Mean)", bold_a)
    worksheet_session_a.write(0, 5, "Peaks Height (Mean)", bold_a)
    worksheet_session_a.write(0, 6, "Peaks Lowpoints (Mean)", bold_a)
    worksheet_session_a.write(0, 7, "Peaks Session (Mean)", bold_a)
    worksheet_session_a.write(0, 8, "Peaks Interval (Mean)", bold_a)
    worksheet_session_a.write(0, 9, "Session Success", bold_a)
    worksheet_session_a.write(0, 10, "Year", bold_a)
    worksheet_session_a.write(0, 11, "Week", bold_a)

    workbook_session_s = xlsxwriter.Workbook(xls_session_s_name)
    worksheet_session_s = workbook_session_s.add_worksheet(name="Main")
    bold_a = workbook_session_s.add_format({'bold': True})
    worksheet_session_s.write(0, 0, "User", bold_a)
    worksheet_session_s.write(0, 1, "Session Count", bold_a)
    worksheet_session_s.write(0, 2, "Duration (Sum)", bold_a)
    worksheet_session_s.write(0, 3, "Duration (Mean)", bold_a)
    worksheet_session_s.write(0, 4, "Peaks (Mean)", bold_a)
    worksheet_session_s.write(0, 5, "Peaks Height (Mean)", bold_a)
    worksheet_session_s.write(0, 6, "Peaks Lowpoints (Mean)", bold_a)
    worksheet_session_s.write(0, 7, "Peaks Session (Mean)", bold_a)
    worksheet_session_s.write(0, 8, "Peaks Interval (Mean)", bold_a)
    worksheet_session_s.write(0, 9, "Session Success", bold_a)
    worksheet_session_s.write(0, 10, "Year", bold_a)
    worksheet_session_s.write(0, 11, "Week", bold_a)



    p = multiprocessing.Pool(processes=multiprocessing.cpu_count() - 1)
    async_result = p.map(worker, files)
    p.close()
    p.join()

    x = 10
    """Write header."""
    bold = workbook.add_format({'bold': True})
    worksheet.write(x - 1, 0, "Chart", bold)
    worksheet.write(x - 1, 1, "User", bold)
    worksheet.write(x - 1, 2, "ID", bold)
    worksheet.write(x - 1, 14, "Date", bold)

    worksheet.write(x - 1, 3, "Year", bold)
    worksheet.write(x - 1, 4, "Filename", bold)
    worksheet.write(x - 1, 5, "Duration", bold)
    worksheet.write(x - 1, 6, "Peaks", bold)
    worksheet.write(x - 1, 7, "Peaks Mean", bold)
    worksheet.write(x - 1, 8, "Lowpoints Mean", bold)
    worksheet.write(x - 1, 9, "SessionMean_<.40", bold)
    worksheet.write(x - 1, 10, "IntervalMean", bold)

    worksheet.write(x - 1, 11, "Peak Widths", bold)
    worksheet.write(x - 1, 12, "PeakHR Session", bold)
    worksheet.write(x - 1, 13, "Success/Failure (1/0)", bold)
    worksheet.write(x - 1, 15, "Peak Timing", bold)


    """Mean stuff"""
    peaks_lowpoints_mean = np.array([y for t in async_result if t is not None for y in t["peaks_lowpoints"]]).mean()
    peaks_values_mean = np.array([y for t in async_result if t is not None for y in t["peaks_values"]]).mean()
    peaks_session_mean = np.array([y for t in async_result if t is not None for y in t["peaks_session"]]).mean()
    peaks_interval_data_mean = np.array([y for t in async_result if t is not None for y in t["peaks_interval_data"]]).mean()

    not_valid = []
    hash_check = {}

    async_result = [x for x in async_result if x is not None]  # Strip none elements

    if sort:
        async_result.sort(key=lambda d: d["date"])

    weekly_stats = {}
    session_data = {}

    for i, res in enumerate(async_result):
        user = res["user"]

        if res is None:
            # Failures...
            continue

        valid = True
        is_duplicate = False
        reasons = []
        duplicate_hash = '|'.join([str(x) for x in res["peak_times"]])  # Just concatinate all peak times
        duplicate_hash += str(res["duration"])
        duplicate_hash += str(res["n_peaks"])
        duplicate_hash += '|'.join([str(x) for x in res["peak_widths"]])  # Just concatinate all peak times

        if duplicate_hash in hash_check:
            valid = False
            reasons.append("DUPLICATE")
            is_duplicate = True
        else:
            hash_check[duplicate_hash] = True

        """Constraints go here.."""
        """if res["n_peaks"] == CONSTRAINT_N_PEAKS:

            if "LOW" in res["peak_types"][1:]:
                valid = False
                reasons.append("LOW_IN_N_PEAKS[1:]")
                print(res["file_name"])
        """

        """if res["n_peaks"] >= CONSTRAINT_N_PEAKS:
            # This is incomplete...
            # Enough peaks
            # Count highs from the back
            high_peaks = res["peak_types"].count("HIGH")
            end_highs = 0
            for e in reversed(res["peak_types"]):
                if e == "LOW":

                    break
                end_highs += 1

            if high_peaks < 3:
                valid = False
                reasons.append("N_HIGH_PEAKS")
        """
        if res["n_peaks"] < CONSTRAINT_N_PEAKS:
            valid = False
            reasons.append("N_PEAKS")


            """n_high = len([x for x in res["peak_types"] if x == "HIGH"])

            if n_high < 3:
                valid = False
                reasons.append("N_PEAKS_HIGH")

            else:
                # To few peaks in general
                valid = False
                reasons.append("N_PEAKS")"""

        if res["duration"] < CONSTRAINT_DURATION:
            valid = False
            reasons.append("DURATION")

        if valid:
            """Move from plots to valid directory"""
            os.rename(
                os.path.join(PLOT_DESTINATION, res["file_name"] + ".png"),
                os.path.join(PLOT_VALID_DESTINATION, res["file_name"] + ".png")
            )

            insert_into_workbook(workbook, worksheet, x, res, valid=True, images=images)
        else:
            not_valid.append(0)
            """Move from plots to invalid directory"""
            if is_duplicate:
                os.remove(os.path.join(PLOT_DESTINATION, res["file_name"] + ".png"))
                x -= 1
            else:
                os.rename(
                    os.path.join(PLOT_DESTINATION, res["file_name"] + ".png"),
                    os.path.join(PLOT_INVALID_DESTINATION, res["file_name"] + ".png")
                )

                insert_into_workbook(workbook, worksheet, x, res, valid=False, images=images, reasons=reasons)
        x += 1

        """
        worksheet.write(x, 1, res["user"])
        worksheet.write(x, 2, str(res["id"]))
        worksheet.write(x, 3, str(res["year"]))
        worksheet.write(x, 4, res["file_name"])
        worksheet.write(x, 5, (str(res["duration"] / 60)))
        worksheet.write(x, 6, (str(res["n_peaks"])))
        worksheet.write(x, 7, str(res["peaks_values"].mean()).replace("nan", "-"))
        worksheet.write(x, 8, str(res["peaks_lowpoints"].mean()).replace("nan", "-"))
        worksheet.write(x, 9, str(res["peaks_session"].mean()).replace("nan", "-"))
        worksheet.write(x, 10, str(res["peaks_interval_data"].mean()).replace("nan", "-"))
        worksheet.write(x, 11, str(res["peak_widths"].mean()).replace("nan", "-"))
        worksheet.write(x, 12, str(res["max_hr"]).replace("nan", "-"))
        worksheet.write(x, 13, int(valid))
        worksheet.write_datetime(x, 14, res["date"], date_format)"""

        if user not in session_data:
            session_data[user] = []

        session_data[user].append([
            user,
            res["id"],
            res["date"], # if 2 then dateformat
            str(res["n_peaks"]),
            str(res["peaks_values"].mean()).replace("nan", "-"),
            str(res["peaks_lowpoints"].mean()).replace("nan", "-"),
            str(res["peaks_session"].mean()).replace("nan", "-"),
            str(res["peaks_interval_data"].mean()).replace("nan", "-"),
            str(res["peak_widths"].mean()).replace("nan", "-"),
            int(valid),
            res["year"]
        ])


        # Weekly statistics TODO
        dt = res["date"]
        week_number = dt.isocalendar()[1]

        if user not in weekly_stats:

            weekly_stats[user] = {}
        if dt.year not in weekly_stats[user]:
            weekly_stats[user][dt.year] = {}
        if week_number not in weekly_stats[user][dt.year]:
            weekly_stats[user][dt.year][week_number] = dict(
                peaks=[],
                duration=[],
                peaks_mean=[],
                peaks_lowpoints=[],
                peaks_session=[],
                peaks_interval_data=[],
                max_hr=[],
                success=[]

            )

        weekly_stats[user][dt.year][week_number]["peaks"].append(res["n_peaks"])
        weekly_stats[user][dt.year][week_number]["duration"].append(res["duration"])
        weekly_stats[user][dt.year][week_number]["peaks_mean"].append(res["peaks_values"].mean())
        weekly_stats[user][dt.year][week_number]["peaks_lowpoints"].append(res["peaks_lowpoints"].mean())
        weekly_stats[user][dt.year][week_number]["peaks_session"].append(res["peaks_session"].mean())
        weekly_stats[user][dt.year][week_number]["peaks_interval_data"].append(res["peaks_interval_data"].mean())
        weekly_stats[user][dt.year][week_number]["max_hr"].append(res["max_hr"])
        weekly_stats[user][dt.year][week_number]["success"].append(int(valid))

        session_data[user].append([
            user,
            res["id"],
            res["date"], # if 2 then dateformat
            str(res["n_peaks"]),
            str(res["peaks_values"].mean()).replace("nan", "-"),
            str(res["peaks_lowpoints"].mean()).replace("nan", "-"),
            str(res["peaks_session"].mean()).replace("nan", "-"),
            str(res["peaks_interval_data"].mean()).replace("nan", "-"),
            str(res["peak_widths"].mean()).replace("nan", "-"),
            int(valid),
            res["year"]
        ])

    # Write the output mean data.
    global_c = 0
    for user, data in session_data.items():
        for rows in data:
            for l, col in enumerate(rows):

                if l == 2:
                    worksheet_session.write_datetime(global_c, l, col, date_format_s)
                else:
                    worksheet_session.write(global_c, l, col)
            global_c += 1



    # Weekly stats
    bold = workbook.add_format({'bold': True})
    o_s_s = 0
    o_s_a = 0
    for _ in range(2):
        # One iteration for successes + failures, another for successes only.

        for user, year_data in weekly_stats.items():

            year_data = collections.OrderedDict(sorted(year_data.items()))
            for year, week_data in year_data.items():
                week_data = collections.OrderedDict(sorted(week_data.items()))


                worksheet.write(x, 0, str(year), bold)
                worksheet.write(x, 3, "Session Count", bold)
                worksheet.write(x, 4, "Weekly Duration", bold)
                if _ == 1:
                    worksheet.write(x, 1, "Successful Only", bold)
                x += 1
                for week, data in week_data.items():


                    if _ == 1:
                        ws = worksheet_session_s
                        o_s_s += 1
                        c = o_s_s
                        indices = np.where(np.array(data["success"]) == 1)
                    else:
                        ws = worksheet_session_a
                        o_s_a += 1
                        c = o_s_a
                        indices = np.array([k for k in range(len(data["success"]))])

                    worksheet.write(x, 0, str(week), bold)
                    worksheet.write(x, 3, str(len(np.array(data["duration"])[indices])))
                    worksheet.write(x, 4, str((np.array(data["duration"])[indices].sum() / 60)).replace("nan", "-")) # kan du lage en akkumulert duration pr uke også? (tenkte da at det letteste hadde vært å legge inn en kolonne til for hele fila, så står denne så klart tom for enkeltfilene som listes først)
                    worksheet.write(x, 5, str((np.array(data["duration"])[indices].mean() / 60)).replace("nan", "-"))
                    worksheet.write(x, 6, str(np.array(data["peaks"])[indices].mean()).replace("nan", "-"))
                    worksheet.write(x, 7, str(np.array(data["peaks_mean"])[indices].mean()).replace("nan", "-"))
                    worksheet.write(x, 8, str(np.array(data["peaks_lowpoints"])[indices].mean()).replace("nan", "-"))
                    worksheet.write(x, 9, str(np.array(data["peaks_session"])[indices].mean()).replace("nan", "-"))
                    worksheet.write(x, 10, str(np.array(data["peaks_interval_data"])[indices].mean()).replace("nan", "-"))
                    worksheet.write(x, 11, str(np.array(data["max_hr"])[indices].mean()).replace("nan", "-"))
                    worksheet.write(x, 12, str(np.array(data["success"])[indices].mean()).replace("nan", "-"))

                    ws.write(c, 0, user)
                    ws.write(c, 1, "-")
                    ws.write(c, 2, str(len(np.array(data["duration"])[indices])))
                    ws.write(c, 3, str((np.array(data["duration"])[indices].sum() / 60)).replace("nan", "-")) # kan du lage en akkumulert duration pr uke også? (tenkte da at det letteste hadde vært å legge inn en kolonne til for hele fila, så står denne så klart tom for enkeltfilene som listes først)
                    ws.write(c, 4, str((np.array(data["duration"])[indices].mean() / 60)).replace("nan", "-"))
                    ws.write(c, 5, str(np.array(data["peaks"])[indices].mean()).replace("nan", "-"))
                    ws.write(c, 6, str(np.array(data["peaks_mean"])[indices].mean()).replace("nan", "-"))
                    ws.write(c, 7, str(np.array(data["peaks_lowpoints"])[indices].mean()).replace("nan", "-"))
                    ws.write(c, 8, str(np.array(data["peaks_session"])[indices].mean()).replace("nan", "-"))
                    ws.write(c, 9, str(np.array(data["peaks_interval_data"])[indices].mean()).replace("nan", "-"))
                    ws.write(c, 10, str(np.array(data["success"])[indices].mean()).replace("nan", "-"))
                    ws.write(c, 11, str(year))
                    ws.write(c, 12, str(week))

                    x += 1


    """worksheet.write(x - 1, 0, "Chart", bold)
    worksheet.write(x - 1, 1, "User", bold)
    worksheet.write(x - 1, 2, "ID", bold)
    worksheet.write(x - 1, 14, "Date", bold)

    worksheet.write(x - 1, 3, "Year", bold)
    worksheet.write(x - 1, 4, "Filename", bold)
    worksheet.write(x - 1, 5, "Duration", bold)
    worksheet.write(x - 1, 6, "Peaks", bold)
    worksheet.write(x - 1, 7, "Peaks Mean", bold)
    worksheet.write(x - 1, 8, "Lowpoints Mean", bold)
    worksheet.write(x - 1, 9, "SessionMean_<.40", bold)
    worksheet.write(x - 1, 10, "IntervalMean", bold)

    worksheet.write(x - 1, 11, "Peak Widths", bold)
    worksheet.write(x - 1, 12, "MaxHR Session", bold)
    worksheet.write(x - 1, 13, "Success/Failure (1/0)", bold)
    worksheet.write(x - 1, 15, "Peak Timing", bold)"""



    """Print a few statistics."""
    worksheet.write(1, 0, "Number of Sessions:", bold)
    worksheet.write(1, 1, str(len(files)))
    worksheet.write(2, 0, "Number of Successful:", bold)
    worksheet.write(2, 1, str(len(files) - len(not_valid)))
    worksheet.write(3, 0, "Number of Failed:", bold)
    worksheet.write(3, 1, str(len(not_valid)))
    worksheet.write(3, 0, "Number of Failed:", bold)
    worksheet.write(3, 1, str(len(not_valid)))

    worksheet.write(4, 0, "Peaks Mean:", bold)
    worksheet.write(4, 1, str(peaks_values_mean).replace("nan", "-"))

    worksheet.write(5, 0, "Lowpoints Mean:", bold)
    worksheet.write(5, 1, str(peaks_lowpoints_mean).replace("nan", "-"))

    worksheet.write(6, 0, "Session Mean < .40:", bold)
    worksheet.write(6, 1, str(peaks_session_mean).replace("nan", "-"))

    worksheet.write(7, 0, "Interval Mean:", bold)
    worksheet.write(7, 1, str(peaks_interval_data_mean).replace("nan", "-"))
    workbook.close()

    workbook_session.close()
    workbook_session_a.close()
    workbook_session_s.close()
    # print(np.array(stats.average_peak_widths).mean())



if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--hrr_low", default=.60, type=float, help="The lower bounds of HRR evaluation.")
    parser.add_argument("--hrr_high", default=.80, type=float, help="The upper bounds of HRR evaluation.")
    parser.add_argument("--filter_hrr", default=.40, type=float,
                        help="Removes datapoints below the given HRR threshold.")
    parser.add_argument("--peak_w_low", default=40, type=int,
                        help="The lower bound of the peak width requirement. In seconds")
    parser.add_argument("--peak_w_high", default=450, type=int,
                        help="The upper bound of the peak width requirement. In seconds")
    parser.add_argument("--minimum_data_points", default=10, type=int,
                        help="Minimum number of datapoints required to be accepted.")
    parser.add_argument("--constraint_duration", default=900, type=int,
                        help="The minimum duration a test must be. Given in seconds.")
    parser.add_argument("--constraint_peaks", default=5, type=int, help="Minimum number of peaks.")
    parser.add_argument("--peak_height", default=10, type=int, help="")

    parser.add_argument("--mode", default="specific", type=str,
                        help="patient_sort, specific, or split. 'patient_sort', generates a separate excel document for each patient with tests in sorted order..'specific' you must also specity the --patient flag with the corresponding ID. ie --patient 008AM")
    parser.add_argument("--patient", default="288DV", type=str,
                        help="When --mode specific is set, this flag specifies which patient to generate.")
    parser.add_argument("--split", default=10, type=int,
                        help="When --mode split is set, this flag specifies how many files to split the generated data into.")
    parser.add_argument("--unzip", default=True, type=bool,
                        help="When this flag is set, zip files are extracted. Takes a few extra seconds, hence this is optional for all runs.")
    parser.add_argument("--output", default="output", type=str, help="Output location")
    parser.add_argument("--dpi", default=100, type=int,
                        help="Resolution of plots. NB Takes longer time when this increases!")
    parser.add_argument("--save_pdf", default=False, type=bool,
                        help="Save as PDF. Takes longer time when this increases!")
    parser.add_argument("--plot_hr_lines", default=True, type=bool, help="")
    parser.add_argument("--plot_width_lines", default=True, type=bool, help="")
    parser.add_argument("--plot_peak_points", default=True, type=bool, help="")
    parser.add_argument("--labeling_mode", default=True, type=bool, help="")
    parser.add_argument("--labeling_count", default=121, type=int, help="")

    args = parser.parse_args()

    CONSTRAINT_DURATION = args.constraint_duration
    CONSTRAINT_N_PEAKS = args.constraint_peaks

    LABELING_MODE = args.labeling_mode
    LABELING_N = args.labeling_count

    """NB. We need to recompile this. remember it is also defined in top of file!"""
    OUTPUT_LOCATION = os.path.join(dir_path, args.output)
    HRM_DESTINATION = os.path.join(OUTPUT_LOCATION, "HRM")
    EXCEL_DEST = os.path.join(OUTPUT_LOCATION, "EXCEL")
    PLOT_DESTINATION = os.path.join(OUTPUT_LOCATION, "PLOTS")
    PLOT_VALID_DESTINATION = os.path.join(PLOT_DESTINATION, "valid")
    PLOT_INVALID_DESTINATION = os.path.join(PLOT_DESTINATION, "invalid")
    LABELING_DESTINATION = os.path.join(OUTPUT_LOCATION, "LABELING")
    os.makedirs(OUTPUT_LOCATION, exist_ok=True)
    os.makedirs(HRM_DESTINATION, exist_ok=True)
    os.makedirs(EXCEL_DEST, exist_ok=True)
    os.makedirs(PLOT_DESTINATION, exist_ok=True)

    shutil.rmtree(PLOT_INVALID_DESTINATION, ignore_errors=True)
    shutil.rmtree(PLOT_VALID_DESTINATION, ignore_errors=True)
    shutil.rmtree(LABELING_DESTINATION, ignore_errors=True)
    os.makedirs(PLOT_INVALID_DESTINATION, exist_ok=True)
    os.makedirs(PLOT_VALID_DESTINATION, exist_ok=True)
    os.makedirs(LABELING_DESTINATION, exist_ok=True)
    """END of mess."""

    util.HRR_PARAM_LOW = args.hrr_low
    util.HRR_PARAM_HIGH = args.hrr_high
    util.END_HRR_STRIP = args.filter_hrr
    util.PEAK_WIDTH_LOW = args.peak_w_low
    util.PEAK_WIDTH_HIGH = args.peak_w_high
    util.TEST_MIN_POINTS = args.minimum_data_points
    util.DPI = args.dpi
    util.SAVE_PDF = args.save_pdf
    util.PLOT_HR_LINES = args.plot_hr_lines
    util.PLOT_PEAK_POINTS = args.plot_peak_points
    util.PLOT_WIDTH_LINES = args.plot_width_lines
    util.PEAK_HEIGHT = args.peak_height
    mode = args.mode
    specific = args.patient
    do_split_n = args.split

    """Unzip and organize files."""
    unzipper(ignore=not args.unzip)

    if mode == "patient_sort":
        chunks = []
        unsorted_files = list(glob.glob(os.path.join(HRM_DESTINATION, "*.hrm")))
        data_per_patient = {}
        for uf in unsorted_files:
            file_name = os.path.splitext(os.path.basename(uf))[0]
            file_name_split = file_name.split("-")
            patient_id = file_name_split[0]
            hrm_id = file_name_split[2]

            if patient_id not in data_per_patient:
                data_per_patient[patient_id] = []

            data_per_patient[patient_id].append([patient_id, hrm_id, uf])

        for k, v in data_per_patient.items():
            data_per_patient[k].sort(key=lambda d: d[1])
            chunks.append([d[2] for d in data_per_patient[k]])

        for chunk in chunks:
            file_name = os.path.splitext(os.path.basename(chunk[0]))[0]
            file_name_split = file_name.split("-")
            pid = file_name_split[0]
            parse(n_start=0, n_files=-1, images=True, override=chunk, prefix=pid)
    elif mode == "specific":
        files = [x for x in list(glob.glob(os.path.join(HRM_DESTINATION, "*.hrm"))) if specific in x]
        parse(n_start=0, n_files=-1, images=True, override=files)
    elif mode == "split":

        n_files = len(list(glob.glob(os.path.join(HRM_DESTINATION, "*.hrm"))))
        per_chunk = int(n_files / do_split_n)

        for x in range(do_split_n):
            start = per_chunk * x
            end = (per_chunk * x) + per_chunk
            parse(n_start=start, n_files=end)

        remaining = n_files - end
        parse(n_start=end, n_files=remaining)
    else:
        parse(n_start=0, n_files=-1, images=False, sort=True)

    if LABELING_MODE:

        INVALID_P = .40
        VALID_P = 1 - INVALID_P
        random.seed("12310")

        invalid_files = os.listdir(PLOT_INVALID_DESTINATION)
        valid_files = os.listdir(PLOT_VALID_DESTINATION)

        if len(invalid_files) >= INVALID_P * LABELING_N and len(valid_files) >= VALID_P * LABELING_N:
            selected_valid = [(os.path.join(PLOT_VALID_DESTINATION, random.choice(valid_files)), "v") for x in
                              range(int(VALID_P * LABELING_N))]
            selected_invalid = [(os.path.join(PLOT_INVALID_DESTINATION, random.choice(invalid_files, )), "i") for x in
                                range(int(INVALID_P * LABELING_N))]

            f = open(os.path.join(LABELING_DESTINATION, "__registry.txt"), "w")

            for i, (labeled_path, t), in enumerate(selected_valid + selected_invalid):
                base_name = os.path.basename(labeled_path)
                labeled_name = "_%s_%s" % (t, base_name)
                hidden_name = "%s.png" % i
                xls_location = base_name.split("-")[0] + ".xls"

                unlabeled_path = os.path.join(PLOT_NO_LABELS_DESTINATION, base_name)

                """Copy the unlabeled sample"""
                shutil.copy(unlabeled_path, os.path.join(LABELING_DESTINATION, hidden_name))

                """Copy the labeled sample"""
                shutil.copy(labeled_path, os.path.join(LABELING_DESTINATION, labeled_name))

                f.write(
                    "--%s--\npng: %s\nxls: %s\n\n" % (hidden_name, labeled_name, xls_location)
                )

            f.close()

        else:
            print("Failed to create labelered data becuase of N")

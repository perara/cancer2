
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

dir_path = os.path.dirname(os.path.realpath(__file__))

pd.set_option('display.expand_frame_repr', False)

PATH_MUST_CONTAIN = ["Hög", "HÖG"]
HR_FILE_EXTENSIONS = ('.hrm')


ZIP_DESTINATION = os.path.join(dir_path, "zips")
OUTPUT_LOCATION = os.path.join(dir_path, "output")
HRM_DESTINATION = os.path.join(OUTPUT_LOCATION, "HRM")
EXCEL_DEST = os.path.join(OUTPUT_LOCATION, "EXCEL")
PLOT_DESTINATION = os.path.join(OUTPUT_LOCATION, "PLOTS")
PLOT_VALID_DESTINATION = os.path.join(PLOT_DESTINATION, "valid")
PLOT_INVALID_DESTINATION = os.path.join(PLOT_DESTINATION, "invalid")

CONSTRAINT_DURATION = 15 * 60
CONSTRAINT_N_PEAKS = 5

os.makedirs(OUTPUT_LOCATION, exist_ok=True)
os.makedirs(HRM_DESTINATION, exist_ok=True)
os.makedirs(EXCEL_DEST, exist_ok=True)
os.makedirs(PLOT_DESTINATION, exist_ok=True)


shutil.rmtree(PLOT_INVALID_DESTINATION, ignore_errors=True)
shutil.rmtree(PLOT_VALID_DESTINATION, ignore_errors=True)

os.makedirs(PLOT_INVALID_DESTINATION, exist_ok=True)
os.makedirs(PLOT_VALID_DESTINATION, exist_ok=True)

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
        hrm_files = [x for x in zip_ref.filelist if x.filename.endswith(HR_FILE_EXTENSIONS) and any(s in x.filename for s in PATH_MUST_CONTAIN)]
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
    #print(hrm_filepath)
    file_name = os.path.basename(hrm_filepath)
    file_name = os.path.splitext(file_name)[0]
    splt = file_name.split("-")
    hrm_user = splt[0]
    hrm_year = splt[1]
    hrm_id = splt[2]

    hrm = HRM(hrm_filepath)
    if hrm.validate():
        hrm.plot(save=True, save_location=PLOT_DESTINATION, save_name=file_name + ".png")
        hrm_peaks = hrm.get_peaks()
        hrm_peak_times = hrm.get_peak_times(hrm_peaks)
        hrm_peak_types = [
            "LOW" if hrm.get_hrr(x) < (util.HRR_PARAM_HIGH*100) else "HIGH" for x in hrm_peaks[0]
        ]

        return dict(
            user=hrm_user,
            year=hrm_year,
            id=hrm_id,
            peak_widths=hrm_peaks[2],
            interval_scale=hrm._hr_interval,
            peak_types=hrm_peak_types,
            duration=hrm.get_duration(),
            n_peaks=len(hrm_peaks[1]),
            peak_times=hrm_peak_times,
            file_name=file_name

        )

    return None


def insert_into_workbook(workbook, worksheet, x, res, valid, images):

    if images:
        sheet_name = str(res["id"]) + str(random.getrandbits(15))
        images_sheet = workbook.add_worksheet(name=sheet_name)
        worksheet.write(x, 10, "=HYPERLINK(\"#%s!B1\", \"Chart\")" % sheet_name)
        plot_dest = PLOT_VALID_DESTINATION if valid else PLOT_INVALID_DESTINATION
        images_sheet.write(0, 1, "=HYPERLINK(\"#Main!K%s\", \"GO BACK\")" % (x + 1))
        images_sheet.insert_image("B2", os.path.join(plot_dest, res["file_name"] + ".png"))

    worksheet.write(x, 0, res["user"])
    worksheet.write(x, 1, str(res["id"]))
    worksheet.write(x, 2, str(res["year"]))
    worksheet.write(x, 3, res["file_name"])
    worksheet.write(x, 4, "%s min" % (str(res["duration"] / 60)))
    worksheet.write(x, 5, "%s peaks" % (str(res["n_peaks"])))
    worksheet.write(x, 6, ', '.join([str(x) for x in np.divide(res["peak_times"], 60)]))
    worksheet.write(x, 15, ', '.join([str(x) for x in np.multiply(res["peak_widths"][0], res["interval_scale"])]))

def parse(n_start=0, n_files=-1, images=True, override=None, prefix="output"):
    stats = Stats()

    """Retrieve files"""
    if override is not None:
        files = override
    else:
        files = list(glob.glob(os.path.join(HRM_DESTINATION, "*.hrm")))
        shuffle(files)
        files = files[n_start:n_files]


    workbook = xlsxwriter.Workbook(os.path.join(EXCEL_DEST, '%s-%s-%s.xlsx') % (prefix, n_start, n_files))

    worksheet = workbook.add_worksheet(name="Main")

    p = multiprocessing.Pool(processes=multiprocessing.cpu_count()-1)
    async_result = p.map(worker, files)
    p.close()
    p.join()

    """Write header."""
    bold = workbook.add_format({'bold': True})
    worksheet.write(2, 0, "User", bold)
    worksheet.write(2, 1, "ID", bold)
    worksheet.write(2, 2, "Year", bold)
    worksheet.write(2, 3, "Filename", bold)
    worksheet.write(2, 4, "Duration", bold)
    worksheet.write(2, 5, "Peaks", bold)
    worksheet.write(2, 6, "Peak Timing", bold)

    x = 3
    not_valid = []
    hash_check = {}
    for i, res in enumerate(async_result):
        if res is None:
            # Failures...
            continue

        valid = True
        reasons = []
        duplicate_hash = '|'.join([str(x) for x in res["peak_times"]])  # Just concatinate all peak times
        duplicate_hash += str(res["duration"])
        duplicate_hash += str(res["n_peaks"])
        duplicate_hash += '|'.join([str(x) for x in res["peak_widths"]])  # Just concatinate all peak times

        if duplicate_hash in hash_check:
            valid = False
            reasons.append("DUPLICATE")
        else:
            hash_check[duplicate_hash] = True

        """Constraints go here.."""
        if res["n_peaks"] == CONSTRAINT_N_PEAKS:

            if "LOW" in res["peak_types"][1:]:
                valid = False
                reasons.append("LOW_IN_N_PEAKS[1:]")
                print(res["file_name"])

        if res["n_peaks"] < CONSTRAINT_N_PEAKS:
            valid = False
            reasons.append("N_PEAKS")

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
            x += 1
        else:
            not_valid.append((res, reasons))
            """Move from plots to invalid directory"""
            os.rename(
                os.path.join(PLOT_DESTINATION, res["file_name"] + ".png"),
                os.path.join(PLOT_INVALID_DESTINATION, res["file_name"] + ".png")
            )

    x += 3
    worksheet.write(x-1, 7, "Failure Reason", bold)
    for res in not_valid:
        insert_into_workbook(workbook, worksheet, x, res[0], valid=False, images=images)
        worksheet.write(x, 7, ','.join(res[1]))
        x += 1

    """Print a few statistics."""
    worksheet.write(1, 0, "Number of Tests:", bold)
    worksheet.write(1, 1, str(len(files)))
    worksheet.write(1, 2, "Number of Successful:", bold)
    worksheet.write(1, 3, str(len(files) - len(not_valid)))
    worksheet.write(1, 4, "Number of Failed:", bold)
    worksheet.write(1, 5, str(len(not_valid)))

    workbook.close()

    #print(np.array(stats.average_peak_widths).mean())


if __name__ == "__main__":

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--hrr_low", default=.60, type=float, help="The lower bounds of HRR evaluation.")
    parser.add_argument("--hrr_high", default=.75, type=float, help="The upper bounds of HRR evaluation.")
    parser.add_argument("--filter_hrr", default=.40, type=float, help="Removes datapoints below the given HRR threshold.")
    parser.add_argument("--peak_w_low", default=40, type=int, help="The lower bound of the peak width requirement. In seconds")
    parser.add_argument("--peak_w_high", default=250, type=int, help="The upper bound of the peak width requirement. In seconds")
    parser.add_argument("--minimum_data_points", default=10, type=int, help="Minimum number of datapoints required to be accepted.")
    parser.add_argument("--constraint_duration", default=900, type=int, help="The minimum duration a test must be. Given in seconds.")
    parser.add_argument("--constraint_peaks", default=5, type=int, help="Minimum number of peaks.")

    parser.add_argument("--mode", default="specific", type=str, help="patient_sort, specific, or split. 'patient_sort', generates a separate excel document for each patient with tests in sorted order..'specific' you must also specity the --patient flag with the corresponding ID. ie --patient 008AM")
    parser.add_argument("--patient", default="008MM", type=str, help="When --mode specific is set, this flag specifies which patient to generate.")
    parser.add_argument("--split", default=10, type=int, help="When --mode split is set, this flag specifies how many files to split the generated data into.")
    parser.add_argument("--unzip", default=True, type=bool, help="When this flag is set, zip files are extracted. Takes a few extra seconds, hence this is optional for all runs.")
    parser.add_argument("--output", default="output", type=str, help="Output location")
    parser.add_argument("--dpi", default=100, type=int, help="Resolution of plots. NB Takes longer time when this increases!")
    parser.add_argument("--save_pdf", default=False, type=bool, help="Save as PDF. Takes longer time when this increases!")


    args = parser.parse_args()

    CONSTRAINT_DURATION = args.constraint_duration
    CONSTRAINT_N_PEAKS = args.constraint_peaks

    """NB. We need to recompile this. remember it is also defined in top of file!"""
    OUTPUT_LOCATION = os.path.join(dir_path, args.output)
    HRM_DESTINATION = os.path.join(OUTPUT_LOCATION, "HRM")
    EXCEL_DEST = os.path.join(OUTPUT_LOCATION, "EXCEL")
    PLOT_DESTINATION = os.path.join(OUTPUT_LOCATION, "PLOTS")
    PLOT_VALID_DESTINATION = os.path.join(PLOT_DESTINATION, "valid")
    PLOT_INVALID_DESTINATION = os.path.join(PLOT_DESTINATION, "invalid")
    os.makedirs(OUTPUT_LOCATION, exist_ok=True)
    os.makedirs(HRM_DESTINATION, exist_ok=True)
    os.makedirs(EXCEL_DEST, exist_ok=True)
    os.makedirs(PLOT_DESTINATION, exist_ok=True)
    shutil.rmtree(PLOT_INVALID_DESTINATION, ignore_errors=True)
    shutil.rmtree(PLOT_VALID_DESTINATION, ignore_errors=True)
    os.makedirs(PLOT_INVALID_DESTINATION, exist_ok=True)
    os.makedirs(PLOT_VALID_DESTINATION, exist_ok=True)
    """END of mess."""

    util.HRR_PARAM_LOW = args.hrr_low
    util.HRR_PARAM_HIGH = args.hrr_high
    util.END_HRR_STRIP = args.filter_hrr
    util.PEAK_WIDTH_LOW = args.peak_w_low
    util.PEAK_WIDTH_HIGH = args.peak_w_high
    util.TEST_MIN_POINTS = args.minimum_data_points
    util.DPI = args.dpi
    util.SAVE_PDF = args.save_pdf

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
        parse(n_start=0, n_files=-1)


import shutil
import subprocess
import time
import zipfile

import psutil

from flask import Flask, request, render_template, send_from_directory, send_file, redirect
import os

from flask_autoindex import AutoIndex

dir_path = os.path.dirname(os.path.realpath(__file__))
app = Flask(__name__)
AutoIndex(app, browse_root=os.path.curdir)

@app.route('/index')
def index():
    return send_from_directory(os.path.join(dir_path, "www",), "index.html")

@app.route('/server_load')
def server_load():
    return render_template("files.html", cpu_usage=psutil.cpu_percent())


def zipdir(path, ziph):
    # ziph is zipfile handle
    for root, dirs, files in os.walk(path):
        if "ARCHIVE" in root:
            continue

        for file in files:
            ziph.write(os.path.join(root, file))

@app.route('/clear_archive', methods=["get"])
def clear_archive():
    full_zip_path = os.path.join(dir_path, "output", "ARCHIVE")
    shutil.rmtree(full_zip_path)
    return redirect("/output", code=302)

@app.route('/archive_current', methods=["get"])
def archive_current():

    full_zip_path = os.path.join(dir_path, "output", "ARCHIVE")
    full_archive_path = os.path.join(dir_path, "output")
    os.makedirs(full_zip_path, exist_ok=True)

    full_zipfile_path = os.path.join(full_zip_path, str(int(time.time())) + ".zip")

    zipf = zipfile.ZipFile(full_zipfile_path, 'w', zipfile.ZIP_DEFLATED)
    zipdir(full_archive_path, zipf)
    zipf.close()

    return redirect("/", code=302)

@app.route('/work/submit', methods=["post"])
def submit_work():
    result = request.form
    command = [v for k, v in result.items() if k == "command"][0]
    args = ["--%s=%s" % (k, v) for k, v in result.items() if k != "command"]

    parent_path = os.path.abspath(os.path.join(dir_path, os.pardir))
    script = parent_path + "/" + command

    subprocess.Popen(["python", script] + args)

    return redirect("/", code=302)

if __name__ == '__main__':
    app.run(debug=True, host="0.0.0.0", port=9292)

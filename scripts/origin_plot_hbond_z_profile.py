"""Plot one H-bond-per-water z profile in Origin's embedded Python.

Run this file from Origin's Code Builder or Python console after generating
``hbond_z_profile.csv`` with ``analyze_hbond_z_profile.py``. Update the paths
below when using a different repository location or analysis output directory.
"""

import csv
import traceback
from pathlib import Path

import originpro as op


REPO_ROOT = Path(r"D:/Tools/AIMD-Analysis-Kit")
INPUT_CSV = REPO_ROOT / "_smoke/raw-data/Cu-aimd1/hbond-z-profile-5001-5003/hbond_z_profile.csv"
OUTPUT_DIR = INPUT_CSV.parent / "origin_png"
WORKBOOK_NAME = "HbondZBook"
GRAPH_NAME = "HbondZGraph"
IMAGE_WIDTH = 2400
LINE_WIDTH = 3


def validate_input(path):
    if not path.is_file():
        raise FileNotFoundError("H-bond z-profile CSV not found: %s" % path)
    with path.open("r", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        header = set(reader.fieldnames or [])
        if not next(reader, None):
            raise ValueError("H-bond z-profile CSV has no data rows: %s" % path)
    required = {"delta_z", "Hbond_per_water", "total_Hbond"}
    if not required.issubset(header):
        raise ValueError("CSV must contain columns %s: %s" % (sorted(required), path))


def main():
    validate_input(INPUT_CSV)
    if op.find_book("w", WORKBOOK_NAME) is not None:
        raise RuntimeError("Origin workbook %s already exists; close it or choose another name" % WORKBOOK_NAME)
    if op.find_graph(GRAPH_NAME) is not None:
        raise RuntimeError("Origin graph %s already exists; close it or choose another name" % GRAPH_NAME)

    op.set_show(True)
    wks = op.new_sheet("w", lname=WORKBOOK_NAME)
    book = wks.get_book()
    if book.name != WORKBOOK_NAME:
        raise RuntimeError("Origin assigned workbook %s instead of %s" % (book.name, WORKBOOK_NAME))
    book.lname = "Cu_aimd1_hbond_z_profile"
    wks.from_file(str(INPUT_CSV), keep_DC=False)
    if wks.rows < 1 or wks.cols < 3:
        raise RuntimeError("Origin imported an incomplete z-profile worksheet")

    graph = op.new_graph(template="line", lname=GRAPH_NAME)
    if graph.name != GRAPH_NAME:
        raise RuntimeError("Origin assigned graph %s instead of %s" % (graph.name, GRAPH_NAME))
    graph.show = True
    layer = graph[0]
    curve = layer.add_plot(wks, coly="B", colx="A", type="line")
    curve.color = "Blue"
    curve.layer.SetNumProp("plot%d.line.width" % (curve.index() + 1), LINE_WIDTH)
    layer.axis("x").title = "Surface-relative z (A)"
    layer.axis("y").title = "H-bonds per water"
    layer.rescale()

    graph.lt_exec("layer -s 1; layer.tickstyle=2;")
    graph.lt_exec('layer -s 1; label -p 50 0 -j 1 -n TITLE "Cu/aimd1 donor-centered H-bond profile";')
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    output = OUTPUT_DIR / "Cu-aimd1-donor-hbond_z_profile.png"
    graph.save_fig(str(output), type="png", replace=True, width=IMAGE_WIDTH)
    if not output.is_file() or output.stat().st_size == 0:
        raise RuntimeError("Origin did not create %s" % output)
    print("Created %s from %s | rows=%d | %s" % (graph.name, INPUT_CSV, wks.rows, output))


if __name__ == "__main__":
    try:
        main()
    except Exception:
        print("FAILED while plotting the H-bond z profile")
        traceback.print_exc()
        raise

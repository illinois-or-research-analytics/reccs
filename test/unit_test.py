import pytest  
from check_outputs_test import run
import subprocess
import os
from pathlib import Path


NETWORK = Path("../data/cen.tsv").resolve()
CLUSTERING = Path("../data/cen_0.01_cm.tsv").resolve()

EVAL_DIR = Path("eval/eval_dir")
STATS_OUTPUT = EVAL_DIR / "end.csv"
REF_OUTPUT = EVAL_DIR / "ref.csv"
REF_DEGREE_SEQUENCES = EVAL_DIR / "ref_degree_sequences.json"
END_DEGREE_SEQUENCE = EVAL_DIR / "end_degree_sequences.json"
SBM_OUTPUT = EVAL_DIR / "sbm.csv"
SBM_DEGREE_SEQUENCE = EVAL_DIR / "sbm_degree_sequences.json"

def run_reccs_pipeline():
    build_dir = Path("../build").resolve()
    original_dir = Path.cwd()
    EVAL_DIR.mkdir(parents=True, exist_ok=True)

    try:
        os.chdir(build_dir)
        subprocess.run(["cmake", ".."], check=True)
        subprocess.run(["make"], check=True)

        with open("test.txt", "w") as out, open("err.txt", "w") as err:
            subprocess.run(["./reccs", str(NETWORK), "-c", str(CLUSTERING), "-o", "output.tsv", "-v"],
                           stdout=out, stderr=err, check=True)

        temp_dirs = sorted(Path(build_dir).glob("temp*/"), key=os.path.getmtime, reverse=True)
        if not temp_dirs:
            raise FileNotFoundError("No temp* directories found")
        clustered_sbm_output = temp_dirs[0] / "clustered_sbm/syn_sbm.tsv"

    finally:
        os.chdir(original_dir)

    subprocess.run(["python3", "extlib/stats.py", "-i", str(build_dir / "output.tsv"),
                    "-e", str(CLUSTERING), "-o", str(STATS_OUTPUT)], check=True)

    subprocess.run(["python3", "extlib/stats.py", "-i", str(clustered_sbm_output),
                    "-e", str(CLUSTERING), "-o", str(SBM_OUTPUT)], check=True)

    subprocess.run(["python3", "extlib/stats.py", "-i", str(NETWORK),
                    "-e", str(CLUSTERING), "-o", str(REF_OUTPUT)], check=True)


def test_RECCS():
    
    run_reccs_pipeline()
    
    assert run(
        stats_dir=STATS_OUTPUT,
        degseq_dir=END_DEGREE_SEQUENCE,
        reference_stats_dir=REF_OUTPUT,
        reference_degseq_dir=REF_DEGREE_SEQUENCES,
        sbm_stats_dir=SBM_OUTPUT,
        sbm_degseq_dir=SBM_DEGREE_SEQUENCE
    ) is True

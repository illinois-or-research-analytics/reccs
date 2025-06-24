import pytest  
import subprocess
import os
from pathlib import Path
import traceback

def run_eval_pipeline():
    try:
        result = subprocess.run(
            ["bash", "eval_pipeline_mini.sh"], 
            capture_output=True, 
            text=True, 
            check=False  # we don't raise exception on non-zero exit
        )

        stdout = result.stdout
        stderr = result.stderr

        print("=== Script Output ===")
        print(stdout)
        print("=== Script Errors ===")
        print(stderr)

        # Check for "FAIL" in either stdout or stderr
        if "FAIL" in stdout or "FAIL" in stderr:
            print("❌ FAIL detected in output.")
            return False
        else:
            print("✅ No FAIL detected.")
            return True

    except Exception as e:
        print(f"🚨 Error running script: {e}")
        return False


def test_RECCS():
        
    assert run_eval_pipeline() is True

import glob
import re
import time
import subprocess

filenames = glob.glob("./S151RfamDataset/*")
MAX_SPAN = 100
TURN = 3

for filename in filenames:
    with open(filename, "r") as f:
        sequence = "".join([x.strip().split(" ")[1] for x in f.readlines()])

        cmd1 = f'./rintp {sequence} {min(MAX_SPAN, len(sequence))} CentroidFold'
        structure = subprocess.check_output(cmd1.split(" ")).decode().strip()

        cmd2 = f'./rintp {sequence} {structure} {min(MAX_SPAN, len(sequence))} RintPwithDFT'
        time_start = time.time()
        subprocess.run(cmd2.split(" "), stdout=subprocess.DEVNULL)
        time_end = time.time()
        print(f"{filename} time(ms)= {int((time_end - time_start)*1000)} length= {len(sequence)}")

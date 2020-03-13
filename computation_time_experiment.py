import glob
import re
import time
import subprocess

filenames = sorted(glob.glob("./S151RfamDataset/*"))
MAX_SPAN = [50, 100, 200, 10000000]
TURN = 3

result = []
count = 0
for filename in filenames:
    with open(filename, "r") as f:
        sequence = "".join([x.strip().split(" ")[1] for x in f.readlines()])
        if re.fullmatch(r"^[ACGU]+$", sequence) is None: continue

        for W in MAX_SPAN:
            actual_W = min(W, len(sequence))
            cmd1 = f'./rintp {sequence} {actual_W} CentroidFold'
            structure = subprocess.check_output(cmd1.split(" ")).decode().strip()

            cmd2 = f'./rintp {sequence} {structure} {actual_W} RintPwithDFT'
            time_start = time.time()
            subprocess.run(cmd2.split(" "), stdout=subprocess.DEVNULL)
            time_end = time.time()
            elapsed_time = int((time_end - time_start)*1000)
            print(f"{count} {filename} time(ms)= {elapsed_time} length= {len(sequence)} W={actual_W}")
            result.append((len(sequence), actual_W, elapsed_time))
            if actual_W == len(sequence): break

with open("computation_time_result.csv", "w", encoding='utf-8') as f:
    f.write("sequence_length,max_span,elapsed_time(ms)\n")
    for x in sorted(result):
        f.write(f"{x[0]},{x[1]},{x[2]}\n")


import glob
import re
import time
import subprocess

filenames = glob.glob("./S151RfamDataset/*")
MAX_SPAN = 100
TURN = 3

for filename in filenames:
    with open(filename, "r") as f:
        data = f.readlines()
        sequence = ""
        bp = []
        for d in [x.strip().split(" ") for x in data]:
            sequence += d[1]

            # if d[2] == "0":
            #     structure += "."
            # elif int(d[0]) < int(d[2]):
            #     if TURN < (int(d[2]) - int(d[0])) and (int(d[2]) - int(d[0])) <= MAX_SPAN:
            #         bp.append((d[1], len(structure)))
            #         structure += "("
            #     else:
            #         structure += "."
            # else:
            #     if TURN < (int(d[0]) - int(d[2])) and (int(d[0]) - int(d[2])) <= MAX_SPAN:
            #         if (bp[-1][0] + d[1]) not in {"AU","CG","GC","GU","UA","UG"}:
            #             structure = structure[:int(bp[-1][1])] + "." + structure[int(bp[-1][1])+1:] + "."
            #         else:
            #             structure += ")"
            #         bp.pop()
            #     else:
            #         structure += "."

        cmd1 = f'./rintp {sequence} {min(MAX_SPAN, len(sequence))} CentroidFold'
        structure = subprocess.check_output(cmd1.split(" ")).decode().strip()

        cmd2 = f'./rintp {sequence} {structure} {min(MAX_SPAN, len(sequence))} RintPwithDFT'
        time_start = time.time()
        subprocess.run(cmd2.split(" "), stdout=subprocess.DEVNULL)
        time_end = time.time()
        print(f"{filename} time(ms)= {int((time_end - time_start)*1000)} length= {len(sequence)}")



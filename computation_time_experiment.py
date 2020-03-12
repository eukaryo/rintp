import glob
import time
import subprocess

filenames = glob.glob("./S151RfamDataset/*")
assert len(filenames) == 151
MAX_SPAN = 100

for filename in filenames:
    with open(filename, "r") as f:
        data = f.readlines()
        sequence = ""
        structure = ""
        bp = []
        for d in [x.strip().split(" ") for x in data]:
            sequence += d[1]
            if d[2] == "0":
                structure += "."
            elif int(d[0]) < int(d[2]):
                bp.append((d[1], len(structure)))
                structure += "("
            else:
                if (bp[-1][0] + d[1]) not in {"AU","CG","GC","GU","UA","UG"} or int(bp[-1][1]) + MAX_SPAN <= len(structure):
                    structure = structure[:int(bp[-1][1])] + "." + structure[int(bp[-1][1])+1:] + "."
                else:
                    structure += ")"
                bp.pop()

        cmd = f'./rintp {sequence} {structure} {min(100,len(sequence))} RintPwithDFT'
        time_start = time.time()
        subprocess.run(cmd.split(" "), stdout=subprocess.DEVNULL)
        time_end = time.time()
        print(f"{filename} {int((time_end - time_start)*1000)} ms")

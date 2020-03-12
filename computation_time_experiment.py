import glob
import time
import subprocess

filenames = glob.glob("./S151RfamDataset/*")
assert len(filenames) == 151

for filename in filenames:
    with open(filename, "r") as f:
        data = f.readlines()
        sequence = ""
        structure = ""
        for d in [x.strip().split(" ") for x in data]:
            sequence += d[1]
            # if d[2] == "0":
            #     structure += "."
            # elif int(d[0]) < int(d[2]):
            #     structure += "("
            # else:
            #     structure += ")"
            structure += "."

        cmd = f'./rintp {sequence} {structure} {min(100,len(sequence))} RintPwithDFT'
        time_start = time.time()
        subprocess.run(cmd.split(" "))
        time_end = time.time()
        print(f"{filename} {int((time_end - time_start)*1000)} ms")


        

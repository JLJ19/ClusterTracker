import os
import numpy as np
import pandas as pd
import openpyxl
import xlwt

nRuns = 1
for i in range(nRuns):
    Command = 'python scripts/run_mot_challenge.py --BENCHMARK MOT20 --SPLIT_TO_EVAL train --TRACKERS_TO_EVAL CTv0-Run'+str(i+1) +' --METRICS HOTA CLEAR --USE_PARALLEL True --NUM_PARALLEL_CORES 4'
    os.system(Command)

Results = []
for i in range(nRuns):
    Results.append(np.loadtxt('data/trackers/mot_challenge/MOT20-train/CTv0-Run'+str(i+1)+'/pedestrian_summary.txt', skiprows=1))
Results = np.array(Results)
df = pd.DataFrame(Results)
filepath = 'MOT20_CTv0_TrainingSummary.xlsx'
df.to_excel(filepath, index=False)

# You can run the following one-liners to get results on a specific tracker from a specific benchmark
# CTv0
#python scripts/run_mot_challenge.py --BENCHMARK MOT20 --SPLIT_TO_EVAL train --TRACKERS_TO_EVAL CTv0 --METRICS HOTA CLEAR --USE_PARALLEL True --NUM_PARALLEL_CORES 6
#python scripts/run_mot_challenge.py --BENCHMARK HT21 --SPLIT_TO_EVAL train --TRACKERS_TO_EVAL CTv0 --METRICS HOTA CLEAR --USE_PARALLEL True --NUM_PARALLEL_CORES 6

# IOU
#python scripts/run_mot_challenge.py --BENCHMARK MOT20 --SPLIT_TO_EVAL train --TRACKERS_TO_EVAL IOU --METRICS HOTA CLEAR --USE_PARALLEL True --NUM_PARALLEL_CORES 6
#python scripts/run_mot_challenge.py --BENCHMARK HT21 --SPLIT_TO_EVAL train --TRACKERS_TO_EVAL IOU --METRICS HOTA CLEAR --USE_PARALLEL True --NUM_PARALLEL_CORES 6

#SORT
#python scripts/run_mot_challenge.py --BENCHMARK MOT20 --SPLIT_TO_EVAL train --TRACKERS_TO_EVAL SORT --METRICS HOTA CLEAR --USE_PARALLEL True --NUM_PARALLEL_CORES 6
#python scripts/run_mot_challenge.py --BENCHMARK HT21 --SPLIT_TO_EVAL train --TRACKERS_TO_EVAL SORT --METRICS HOTA CLEAR --USE_PARALLEL True --NUM_PARALLEL_CORES 6

#OC-SORT
#python scripts/run_mot_challenge.py --BENCHMARK MOT20 --SPLIT_TO_EVAL train --TRACKERS_TO_EVAL OCSORT --METRICS HOTA CLEAR --USE_PARALLEL True --NUM_PARALLEL_CORES 6
#python scripts/run_mot_challenge.py --BENCHMARK HT21 --SPLIT_TO_EVAL train --TRACKERS_TO_EVAL OCSORT --METRICS HOTA CLEAR --USE_PARALLEL True --NUM_PARALLEL_CORES 6

# ClusterTracker

Description:

This repository contains the code used for ClusterTracker, a multiple object tracking method that is listed as "CTv0" in the MOT20 and HT21 challenges of https://motchallenge.net.  ClusterTracker uses DBSCAN to perform spatial clustering on the object detections from short 10-frame sequences, treats clusters as tracklets, and then merges successive tracklets with high bounding box overlap to form tracks.  ClusterTracker can be run in an online mode with a 10 frame latency or in parallel across multiple cores for offline analysis. ClusterTracker is described in detail in our paper: https://papers.ssrn.com/sol3/papers.cfm?abstract_id=4102945
We hope that this work is useful to you!  We have a few instructions to try to help you get everything working.

Requirements:

Install Python for MATLAB...we used Python 3.8 and MATLAB 2021a (with additional official toolboxes, i.e., Computer Vision, Image Processing, Stats and ML, and Signal Processing).  MATLAB has trouble with conda/virtual environments etc, so you will likely need to install Python 3.8 and then install the additional requirements for the basic functionality of IOU, SORT, and OCSORT - which I believe is just numpy.  You need to then make sure MATLAB is pointing to Python's location...instructions below, "pyenv" in Command Window tells you which python version is active and where it is looking.

https://www.mathworks.com/help/matlab/matlab_external/install-the-matlab-engine-for-python.html

I stripped out most unnecessary code from IOU, SORT, and OCSORT for the purposes of making a barebones example.  Please visit the individual repositories for more info on each tracker - original licenses are included in each tracker's folder.

IOU
https://github.com/bochinski/iou-tracker

SORT
https://github.com/abewley/sort

OCSORT
https://github.com/noahcao/OC_SORT

Getting Started

Once you can get MATLAB to see the correct version of python, in MATLAB do the following:
1. Open up the IOU folder - run py.iou_tracker.TestFunc() in the command window
	1(a) If everything is working, then the message "MATLAB sees IOU" will appear in the command window.
2. Open up the SORT folder - run py.run_sort_once.TestFunc() and get a similar result
3. Open up the OCSORT folder - run py.run_ocsort_once.TestFunc() and get a similar result.
4. Steps 1-3 create a py_cache folder that MATLAB needs.  If this step is not done, then MATLAB will give errors pretty much saying it doesn't recognize your py.stuff.  This may need to be done every time to re-start MATLAB.
5. Navigate in MATLAB to the root folder for this project. Type "addpath IOU SORT OCSORT" in the CommandWindow
6. Type "savepath" in CommandWindow so that MATLAB remembers step 5.
7. Run one of the following scripts.  Scripts not mentioned are function files.

Running Scripts:

MOT20_RunAllTrackers.m

This script allows you to provide settings for each tracker and run all of them on on a MOT Challenge Style Dataset.  We provide MOT20 and HT21 in the Datasets folder - no images.  Note, the MOT20 detection files only contain 0 or 1 for confidence.  Typically, MOT methods will do some sort of adaptive filtering based on detection confidence.  In this script, all detection confidence values are set to 1.

MOT20_Timing.m and HT21_Timing.m

These scripts run all three trackers 11 times (changeable) and write the timing results to a text file for each tracker.  Running them 11 times can take a couple hours. Ignore the first run, because sometimes this can be a lot slower than the rest in MATLAB.  We usually just drop the first row from each text file because of this.

MOT20_ParameterSearch_CTv0.m and HT21_ParameterSearch_CTv0.m

This script runs CTv0 on a custom range of combinations of eps and IOU_thresh.  It writes the results of each combination to a new file, e.g., CTv0-Run1.  3 entries of eps and 3 entries of IOU_thresh will result in 9 combinations.  If you have the official track evaluation code from 
https://github.com/JonathonLuiten/TrackEval, you can treat each run as a different tracker and run our EvalMultipleTrackers.py script from the TrackEval root folder to summarize these results.

EvalMultipleTrackers.py

Run from root folder of the official track evaluation code from 
https://github.com/JonathonLuiten/TrackEval
Helps you group all results from CTv0-Run1 to CTv0-RunX.  Also contains the one-line text to evaluate individual trackers.

BlobTracking_MethodComparison.m

This script allows you to re-create the background-subtraction example from the paper.  Rather than using an object detection file, this script performs background subtraction using a MATLAB tool and saves all the foreground detections. All trackers can be run with these foreground detections.  Then a few filters are applied to the results from each tracker to remove tracks of non-objects.


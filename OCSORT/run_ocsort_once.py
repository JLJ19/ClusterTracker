# Matlab wrapper function for running OC_SORT on a MOT20 style detection file
# that is sent to this function as an (n, 10) size np.array - n= total number of detections

import time
import numpy as np
from trackers.ocsort_tracker.ocsort import OCSort

def TestFunc():
    print('MATLAB sees OCSORT')

def run_matlab_wrapper(seq_trks, var_det_thresh, var_max_age, var_min_hits, var_iou_thresh, var_deltat, var_inertia):
    cats = ["Pedestrian"] # Not sure if this does anything, just avoiding errors
    total_time = 0 
    total_frame = 0 
    # OCSORT settings for MOT20 file using all detections
    # var_det_thresh=-0.1, 
    # var_max_age=30, 
    # var_min_hits=3, 
    # var_iou_thresh=0.3,
    # var_deltat=3, 
    # var_inertia=0.2):
    tracker = OCSort(
    det_thresh=float(var_det_thresh), 
    max_age=int(var_max_age), 
    min_hits=int(var_min_hits),
    iou_threshold=float(var_iou_thresh),
    delta_t=int(var_deltat),
    inertia=float(var_inertia))
    

    # seq_trks is the mot detection matrix
    out = []
    min_frame = seq_trks[:,0].min()
    max_frame = seq_trks[:,0].max()
    for frame_ind in range(int(min_frame), int(max_frame)+1):
        dets = seq_trks[np.where(seq_trks[:,0]==frame_ind)][:,2:6]
        cates = np.zeros((dets.shape[0],))
        scores = seq_trks[np.where(seq_trks[:,0]==frame_ind)][:,6]
        dets[:, 2:] += dets[:, :2] # xywh -> xyxy
        assert(dets.shape[0] == cates.shape[0])
        #t0 = time.time()
        t0 = time.perf_counter()
        online_targets = tracker.update_public(dets, cates, scores)
        #t1 = time.time()
        t1 = time.perf_counter()
        total_frame += 1
        total_time += t1 - t0
        trk_num = online_targets.shape[0]
        boxes = online_targets[:, :4]
        ids = online_targets[:, 4]
        frame_counts = online_targets[:, 6]
        sorted_frame_counts = np.argsort(frame_counts)
        frame_counts = frame_counts[sorted_frame_counts]
        cates = online_targets[:, 5]
        cates = cates[sorted_frame_counts].tolist()
        cates = [cats[int(catid)] for catid in cates]
        boxes = boxes[sorted_frame_counts]
        ids = ids[sorted_frame_counts]
        for trk in range(trk_num):
            lag_frame = frame_counts[trk]
            if frame_ind < 2*var_min_hits and lag_frame < 0:
                continue
            out.append(np.array([int(frame_ind+lag_frame), int(ids[trk]),
            boxes[trk][0], boxes[trk][1], 
            boxes[trk][2]-boxes[trk][0],
            boxes[trk][3]-boxes[trk][1], 1,-1,-1,-1]))
    out = np.array(out) # make np.array out of entire results file
    return total_time, out;
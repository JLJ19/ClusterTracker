% This function passes OCSORT arguments to the matlab_wrapper for OCSORT in
% python.
% py.run_ocsort_once.TestFunc(); % run this to debug and see if
 % Matlab is connecting to python
function [track_result, Time] = run_OCSORT(Dets, OCSORT_Settings)
    var_det_thresh          = OCSORT_Settings.det_thresh; 
    var_max_age             = OCSORT_Settings.max_age; % paper did not comment on this value, but 30=default
    var_min_hits            = OCSORT_Settings.min_hits; % paper did not comment on this value, but 3=default
    var_iou_threshold       = OCSORT_Settings.iou_threshold;
    var_delta_t             = OCSORT_Settings.delta_t;
    var_inertia             = OCSORT_Settings.inertia;

    try
        ret = py.run_ocsort_once.run_matlab_wrapper(py.numpy.array(Dets), var_det_thresh, var_max_age, var_min_hits, var_iou_threshold, var_delta_t, var_inertia);
    catch exception
        disp('error while calling the python tracking module: ')
        disp(' ')
        disp(getReport(exception))
    end
    Time = ret{1};
    track_result = double(ret{2}); % [f, ID, TLx, TLy, W, H, Confidence, -1, -1, -1
end

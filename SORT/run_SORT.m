function [track_result, Time] = run_SORT(Dets, SORT_Settings)
max_age         = SORT_Settings.max_age;
min_hits        = SORT_Settings.min_hits;
iou_threshold   = SORT_Settings.iou_threshold;
try
    ret = py.run_sort_once.SORT_MatlabWrapper(py.numpy.array(Dets), max_age, min_hits, iou_threshold);
catch exception
    disp('error while calling the python tracking module: ')
    disp(' ')
    disp(getReport(exception))
end
Time = ret{1};
track_result = double(ret{2}); % [f, ID, TLx, TLy, W, H, Confidence
end



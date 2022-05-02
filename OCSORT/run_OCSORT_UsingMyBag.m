% This function passes OCSORT arguments to the matlab_wrapper for OCSORT in
% python.
function [Bag, Tubes, Time] = run_OCSORT_UsingMyBag(Bag, MinDuration, OCSORT_Settings)
    Dets = [Bag(:,3) -1*ones(length(Bag(:,1)),1) Bag(:,4:7) ones(length(Bag(:,1)),1) -1*ones(length(Bag(:,1)),3)];
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
    track_result(:,7:end) = [];
    % Change format to [Cx, Cy, f, TLx, TLy, W, H, Sz, id]
    track_result = [track_result(:,3:4)+track_result(:,5:6)/2, track_result(:,1), ...
                    track_result(:,3:6),track_result(:,5).*track_result(:,6), track_result(:,2)];

    % This step filters out any tracks that are too short, and then also resets
    % the track IDs to start at small numbers.
    TubeFlag = false;
    IDs = unique(track_result(:,end));
    j=1;
    for i=1:length(IDs)
        temp = track_result(track_result(:,end) == IDs(i), :);
        if length(temp(:,1)) >= MinDuration
            Tubes{j} = cat(2,temp(:,1:8), j*ones(length(temp(:,1)),1));
            j=j+1;
            TubeFlag = true;
        end
    end
    if TubeFlag==true
        Bag = round(cell2mat(Tubes'));
    else
        Tubes = {};
        Bag = [];
    end
end




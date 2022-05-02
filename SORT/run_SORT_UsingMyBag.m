function [Bag, Tubes, Time] = run_SORT_UsingMyBag(Bag, MinDuration, SORT_Settings)
%Bag = [Cx, Cy, f, TLx, TLy W, H, Sz] - input format
% SORT expects detection in MOTC format ...% Dets = [f, -1, TLx, TLy, W, H, Confidence
Dets = [Bag(:,3) -1*ones(length(Bag(:,1)),1) Bag(:,4:7) ones(length(Bag(:,1)),1) -1*ones(length(Bag(:,1)),3)];
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
track_result(:,7:end) = [];
% Change format to [Cx, Cy, f, TLx, TLy, W, H, Sz, id]
track_result = [track_result(:,3:4)+track_result(:,5:6)/2, track_result(:,1), ...
                track_result(:,3:6),track_result(:,5).*track_result(:,6), track_result(:,2)];

% This step filters out any tracks that are too short, and then also resets
% the track IDs to start at small numbers.
IDs = unique(track_result(:,end));
j=1;
for i=1:length(IDs)
    temp = track_result(track_result(:,end) == IDs(i), :);
    if length(temp(:,1)) >= MinDuration
        Tubes{j} = cat(2,temp(:,1:8), j*ones(length(temp(:,1)),1));
        j=j+1;
    end
end
Bag = round(cell2mat(Tubes'));




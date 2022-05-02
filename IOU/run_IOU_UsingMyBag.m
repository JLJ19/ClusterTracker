 function [Bag,Tubes, Time] = run_IOU_UsingMyBag(Bag, MinTrackDuration, IOU_Settings)
baselinedetections = [Bag(:,3) -1*ones(length(Bag(:,1)),1) Bag(:,4:7) ones(length(Bag(:,1)),1)];

sigma_l     = IOU_Settings.sigma_l;
sigma_h     = IOU_Settings.sigma_h;
sigma_iou   = IOU_Settings.sigma_iou;
t_min       = IOU_Settings.t_min;

%% CVPR19 Settings (Pretty much same as MOT20)
% sigma_l = 0;
% sigma_h = 0.1;
% sigma_iou = 0.3;
% t_min = 1;


% HT21 Settings (From Different Paper)
% sigma_l = 0;
% sigma_h = 0.1;
% sigma_iou = 0.3;
% t_min = 2;


% Mask R-CNN (frcnn)
% sigma_l = 0;
% sigma_h = 0.95;
% sigma_iou = 0.6;
% t_min = 7;

% %% R-CNN
% sigma_l = 0;
% sigma_h = 0.7;
% sigma_iou = 0.5;
% t_min = 2;

% %% ACF
% sigma_l = 0;
% sigma_h = 0.3;
% sigma_iou = 0.5;
% t_min = 3;

% %% CompACT
% sigma_l = 0;
% sigma_h = 0.2;
% sigma_iou = 0.4;
% t_min = 2;

% %% EB
% sigma_l = 0;
% sigma_h = 0.8;
% sigma_iou = 0.5;
% t_min = 2;

%% running tracking algorithm
try
    ret = py.iou_tracker.track_iou_matlab_wrapper(py.numpy.array(baselinedetections(:).'), sigma_l, sigma_h, sigma_iou, t_min);
catch exception
    disp('error while calling the python tracking module: ')
    disp(' ')
    disp(getReport(exception))
end
Time = ret{1};
track_result = cell2mat(reshape(ret{2}.cell.', 6, []).'); % BBox, frame, id

% Change format to [Cx, Cy, f, TLx, TLy, W, H, Sz, id]
track_result = [track_result(:,1:2)+track_result(:,3:4)/2, ...
                track_result(:,5), track_result(:,1:4), ...
                track_result(:,3).*track_result(:,4), track_result(:,end)];

% This step filters out any tracks that are too short, and then also resets
% the track IDs to start at small numbers.
IDs = unique(track_result(:,end));
j=1;
for i=1:length(IDs)
    temp = track_result(track_result(:,end) == IDs(i), :);
    if length(temp(:,1)) >= MinTrackDuration
        Tubes{j} = cat(2,temp(:,1:8), j*ones(length(temp(:,1)),1));
        j=j+1;
    end
end
Bag = round(cell2mat(Tubes'));



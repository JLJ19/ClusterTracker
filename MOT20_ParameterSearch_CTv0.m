% Joachim Lohn-Jaramillo
% 4/19/22

% This script uses CTv0 with a range of parameters.  This is used
% to find the the best set of results.  Each run is output to the Results
% folder as CTv0-Run1, CTv0-Run2, etc.  These files can be copied and
% pasted into the official MOT Evaluation code found here: https://github.com/JonathonLuiten/TrackEval
% To evaluate, paste them into the following path from TrackEval...
%           TrackEval-master\data\trackers\mot_challenge\MOT20-train
% Then the included python script can be run from the TrackEval root folder
% to automatically evaluate all runs, and summarize the results into one excel file.
% To figure out what settings go with what run, a text file is written
% called "CTv0_ParameterList" - row 1 goes with CTv0-Run1, etc.


% -------------------------------------------------------------------------
% User settings

clc
clear all
close all

% Get filenames from folder that we will analyze
Directory = 'Datasets/MOT20/train';
Files = dir(Directory);
Files = Files(3:length(Files)) % Get rid of trash entries - Windows version
% -------------------------------------------------------
% -------------------------------------------------------------------------
% Cluster Tracker Settings
% -----------------------
% CT_v0 Settings
CT_Settings.BinWidth = 10; % Could be increased in multiples of 10
CT_Settings.MaxTubeDim = 1000; % Max bbox dimension allowed for an object
CT_Settings.minpts = 2; % DBSCAN parameter

% BEST MOT20 MOTA params, eps=16, MinOverlapRatio = 0.4, MinBins = 2
IOU_Thresh = [0 0.4]; % intersection/union - used to join tubelets
eps = [16 18]; %
MinBins = 2;

DetConf_Thresh = -0.1;
% -------------------------------------------------------------------------
%% Read detections from file 
% Det = [f, id, TLx, TLy, W, H, Conf]
Parameters = [];
for m=1:length(eps)
    for n=1:length(IOU_Thresh)
        Name = sprintf('Output/MOT20-train/CTv0-Run%s/data', num2str((m-1)*length(IOU_Thresh)+n));
        mkdir(Name)
        Parameters = [Parameters; [(m-1)*length(IOU_Thresh)+n, eps(m), IOU_Thresh(n)]];
    end
end
writematrix(Parameters, 'Output/MOT20-train/CTv0_ParameterList.txt') ;
for m=1:length(eps)
    for n=1:length(IOU_Thresh)
        fprintf('Test: %s\n',num2str((m-1)*length(IOU_Thresh)+n))
        for i=1:length(Files)
            % -------------------------------------------------------------------------
            name{i} = Files(i).name;
            disp(name{i})
            InputVidName = sprintf('%s.mp4', name{i}(1:8));
            Det = readmatrix(sprintf('%s/%s/det/det.txt', Directory, name{i}));
            Det = sortrows(Det, 1);
            nFrames = max(Det(:,1));
            
            Det = Det(Det(:,7) >  DetConf_Thresh, :);
            % Change format to [Cx, Cy, f, TLx, TLy, W, H, Sz]
            FV = [Det(:,3:4)+Det(:,5:6)/2, ...
                            Det(:,1), Det(:,3:6), ...
                            Det(:,5).*Det(:,6)];
            % ---------------------------------------------------------------------
            % ---------------------------------------------------------------------
            %% CT-v0
            CT_Settings.eps = eps(m); 
            CT_Settings.TrackletMerge_OverlapRatio = IOU_Thresh(n); % "IOU Thresh" above
            CT_Settings.MinBins = 2;
            [CT_Tracks, StartPt, EndPt, Stationary, Time_CTv0((m-1)*length(IOU_Thresh)+n)] = CT_v0(FV, nFrames, 0, CT_Settings);
            Bag = round(cell2mat(CT_Tracks'));
            % ---------------------------------------------------------------------
            SaveTracksMotStyle(Bag, sprintf('Output/MOT20-train/CTv0-Run%s/data/%s.txt', num2str((m-1)*length(IOU_Thresh)+n),name{i}))
        % -------------------------------------------------------------------------
        end
    end
end

% -------------------------------------------------------------------------
function [] = SaveTracksMotStyle(Bag, Name)
    % ---------------------------------------------------------------------
    % MyBAG =  [Cx Cy f TLx TLy Width Height TubeID]
    % MOT Data Style = <frame>, <id>, <bb_left>, <bb_top>, <bb_width>, <bb_height>, <conf>, <x>, <y>, <z>
    MotBag_Det = [Bag(:,3), Bag(:,end), Bag(:,4:7), ones(size(Bag(:,1))), -1*ones([length(Bag(:,1)) 3])];
    writematrix(MotBag_Det, Name);
    % ---------------------------------------------------------------------
end



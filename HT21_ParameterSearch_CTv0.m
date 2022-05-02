% JLJ
% This script evaluates a tracker with a range of parameters.  This is used
% to find the the best set of results

% -------------------------------------------------------------------------
% User settings

clc
clear all
close all

% Get filenames from folder that we will analyze
Directory = 'data/HT21/train/';
Files = dir(Directory);
Files = Files(3:length(Files)) % Get rid of trash entries - Windows version
% -------------------------------------------------------

MakeVideos = false;
ShowTubeID = true;
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% Cluster Tracker Settings
% -----------------------
% CT_v0 Settings
CT_Settings.BinWidth = 10; % Could be increased in multiples of 10
CT_Settings.BinOverlap = 0.5; % DONT CHANGE - WOULD PROBABLY BREAK THINGS!!!
CT_Settings.MaxTubeDim = 1000; % Max bbox dimension allowed for an object
CT_Settings.minpts = 2; % DBSCAN parameter, DONT CHANGE

% BEST MOT20  params, eps=20, MinOverlapRatio = 0.60, MinTubeLength = 10
IOU_Thresh = [0 0.4]; 
eps = [16 18];
MinBins = 2;
DetConf_Thresh = 59;


% VisualizeDetStats(Files,Directory)
% -------------------------------------------------------------------------
%% Read detections from file 
% Det = [f, id, TLx, TLy, W, H, Conf]
Parameters = [];
for m=1:length(eps)
    for n=1:length(IOU_Thresh)
        Name = sprintf('Results/HT21-train/CTv0-Run%s/data', num2str((m-1)*length(IOU_Thresh)+n));
        mkdir(Name)
        Parameters = [Parameters; [(m-1)*length(IOU_Thresh)+n, eps(m), IOU_Thresh(n)]];
    end
end
writematrix(Parameters, 'Results/HT21-train/CTv0_ParameterList.txt') ;
for m=1:length(eps)
    for n=1:length(IOU_Thresh)
        fprintf('Test: %s\n',num2str((m-1)*length(IOU_Thresh)+n))
        for i=1:length(Files)
            % -------------------------------------------------------------------------
            name{i} = Files(i).name;
            disp(name{i})
%             InputVidName = sprintf('%s.mp4', name{i}(1:8));
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
            % ---------------------------------------------------------------------
            %% CT-v0
            CT_Settings.eps = eps(m); 
            CT_Settings.TrackletMerge_OverlapRatio = IOU_Thresh(n); % "IOU Thresh" above
            CT_Settings.MinBins = MinBins;
            [CT_Tracks, StartPt, EndPt, Stationary, TrainTime_CTv0((m-1)*length(IOU_Thresh)+n)] = CT_v0(FV, nFrames, 0, CT_Settings);
            Bag = round(cell2mat(CT_Tracks'));
            % ---------------------------------------------------------------------
            SaveTracksMotStyle(Bag, sprintf('Results/HT21-train/CTv0-Run%s/data/%s.txt', num2str((m-1)*length(IOU_Thresh)+n),name{i}))
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





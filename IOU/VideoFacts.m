% Analyze gt files

clc
clear all
close all

txt = readcell('seqmaps/JLJ-train.txt')
txt(1) = [];

Facts = [];

for i=1:length(txt)
    Targets = [];
    Data = readmatrix(sprintf('Custom_Labels/train/%s/gt/gt.txt',txt{i}));
    % Data = [f,id, bbox, conf, ..]
    for j=1:max(Data(:,1))
        tempFrame = Data(Data(:,1) == j,:);
        Targets(j) = length(tempFrame(:,1));
    end
    MeanTargetsPerFrame(i) = mean(Targets);
end


clear all

% List of files
fileList = {'TMR_192.mat';'TMR_193.mat';'TMR_194.mat'};

% List of fields to use
flds = who('-file','TMR_192.mat');
flds(ismember(flds,'diss_TMR_192')) = []; % Get rid of Rolf's fucking stupid named structure

% Create an instance of the MatFiles object.
[mfo] = MatFiles(fileList, flds);
mfo = SetAlias(mfo,'APz','NEMOAPZ');

a = mfo.Ay(921598:921602,1)

% Look at it
% methods(mfo)
% properties(mfo)
% a = fieldnames(mfo)
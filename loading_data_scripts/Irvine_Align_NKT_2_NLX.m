%%%% Script to aign the Nihon Kohden Data and Neuralynx data %%%%
% At UC Irvine, there are two recording systems that must be aligned in
% order to use the data. The first is the Nihan Kohden data, which is the
% system that records the data used in any clinical decisions and often has
% more electrodes. Additionally, there are Neuralynx data which are set up
% by the Irvine team, often have single units, more noise, less channels,
% and the photodiode channel. In order to use the Nihan Kohden data, you
% have to do a cross correlation anaylsis to line up the nerual signals in
% time so that the photodiode channel is usable. This script identifies the
% proper way to lag the Neuralynx (NLX) data so that the photodiode is
% usable for either dataset.

% Be careful with the script and any expectations you may have between the
% two datasets. Sometimes the +/- of the two datasets are switched,
% sometimes electrodes are mislabeled, there are breaks in the NLX data,
% etc. Don't run this straight through and take whatever number it returns.
% Check each step to make sure it makes sense and please add to it as you
% see improvements!

% credit to Ludo for writing the skeleton, Brooke orgnaized it and added
% some comments :)


%% add fieldtrip
ftDir = '~/Projects/fieldtrip-20191213';
if exist('ft_defaults.m', 'file') == 0
addpath(ftDir); ft_defaults;
end

%% useful paths %%
addpath(genpath('/Users/bstavel/Projects/knight_server/remote/lbellier/DataWorkspace/_tools/NLX2mat/'));
addpath('/Users/bstavel/Projects/knight_server/remote/lbellier/DataWorkspace/_tools/');
addpath('/Users/bstavel/Projects/knight_server/remote/lbellier/DataWorkspace/_projects/auditoryLocalizer/');

%% paths to NLX and NKT data
pathNLX = '/Users/bstavel/Projects/pacman/ieeg_data/IR101/raw_data/su_data/';
pathNKT = '/Users/bstavel/Projects/pacman/ieeg_data/IR101/raw_data/nihan_kodan_data/EA0291BL_15-Sep-2021_100951-121021.mat';

%% load NLX data
cfg = [];
cfg.dataset = pathNLX;
dataSU = ft_preprocessing(cfg);

%% load NKT data
load(pathNKT)

%% find the correct segment of Neuralynx data
% there is some strange way that nauralynx records data-- it collects a
% certain set of samples per unit of time that it saves as rows instead of
% columns. There are sometimes breaks in this system where a timepoint is
% missing a certain sample. Hopefully, this doesn't happen in the middle of
% the task. If it does, talk with Ludo. The code identifies the breaks and
% then redfines the trial appropriately. Missing some of the logic on how
% to know which part of the data has your task data since there is no
% photdiode with the nlx data

% load one of the nlx channels and look for break points
ncs = read_neuralynx_ncs(sprintf('%s/LAM1.ncs', pathNLX)); % pick any su channel in your pathNLX data path
idx = find(ncs.NumValidSamp<512);
nValSamp = ncs.NumValidSamp(idx);
nTrials = length(idx);
trl = zeros(nTrials, 3);
trl(:, 2) = ((idx-1)*512 + nValSamp)';
trl(1, 1) = 1;
trl(2:end, 1) = (idx(1:end-1)*512+1)';

% select the segment that contains your task
trl % see how many breaks there are in your segment of data
% idxSUsegment = 1
idxSUsegment = 2
trl = trl(idxSUsegment, :);

% redefine the trial based on chosen chunk
cfg = [];
cfg.trl = trl;
dataSU = ft_redefinetrial(cfg, dataSU);
dataSU.time{1} = dataSU.time{1}(1:size(dataSU.trial{1},2));
dataSU.sampleinfo = [1 length(dataSU.time{1})];
pos1 = [1 -516 1144 1200]; % to resize ft_databrowser

% run to check if the new trial time makes sense
cfg = [];
cfg.viewmode = 'vertical';
cfg.preproc.demean = 'yes';
cfg.blocksize = 10;
cfg.position = pos1;
cfg.ylim = [-1 1].*400;
ft_databrowser(cfg, dataSU);

% if you are happy with this, save the cut point
nlx_start_point = trl(1);

%% create indices for reordering misaslinged electrodes %%
% this bit of code _helps_ you to create an apples to apples comparison
% between the nlx and nkt data by getting indices that match up the electrode
% names across the system. *Importantly*, sometimes the electrodes are
% misnamed on the nlx system... in this case you will have to change the
% correlation approach to look not just across time, but across electrodes


% make dictionary between the two electrodes
NKTlabels = sort(data.header.label);
NLXlabels = sort(dataSU.label);
NLXtoNKT = {'LAM', 'AMG'; 'LHH', 'HH'; 'LIN', 'AIN'; 'LOF', 'OF'; 'LTH', 'TH'; 'LFO', 'PTR'};

% labels to grep on
nlx_grep_labels = cellfun(@(x) x(regexp(x, '\D')), NLXlabels, 'UniformOutput', false)
nkt_grep_labels = cellfun(@(x) x(regexp(x, '\D')), NKTlabels, 'UniformOutput', false)

% get indices of dictionary for the grepped labels
[~, idxElecs1] = ismember(nlx_grep_labels, NLXtoNKT(:, 1)); % NLX
[~, idxElecs2] = ismember(nkt_grep_labels, NLXtoNKT(:, 2)); % NKT

% order and sort the nlx elecs %
[a b] = sort(idxElecs1);
ordered_nlx_elecs = NLXlabels(b);
usable_nlx_elecs = ordered_nlx_elecs(a > 0);

% order and sort the nkt elecs %
[a b] = sort(idxElecs2);
ordered_nkt_elecs = NKTlabels(b);
usable_nkt_elecs = ordered_nkt_elecs(a > 0);

% load NKT data
dataTMP = load(pathNKT);
hdr = dataTMP.data.header;
L_nkt = hdr.nSamples;
taskTimeLoc = [round(L_nkt/2) L_nkt]; % use the second half of the file, from 11:09:51 to 12:10:21, task being from 11:29 to 11:54:18
elecsToRemove = {'E', 'DC01', 'DC02', 'DC03', 'DC04', 'Mark1', 'Mark2', 'X32', 'X33', 'REG', 'EKG', 'Events'};
[~, idxBad] = ismember(hdr.label, elecsToRemove);
[~, idxBad] = ismember(elecsToRemove, hdr.label);
idxOk = setdiff(1:length(hdr.label), idxBad);
dataNKT = ft_array2data(dataTMP.data.data.data(idxOk, taskTimeLoc(1):taskTimeLoc(2))', hdr.Fs, hdr.label(idxOk));
fs = dataNKT.fsample;
clear dataTMP

% cfg = [];
% cfg.trl = [395834 1475360 0];
% dataNKT = ft_redefinetrial(cfg, dataNKT);

cfg = [];
cfg.dataset = pathNLX;
cfg.trl = [1 nlx_start_point 0];
dataNLX = ft_preprocessing(cfg);
fsSU = dataNLX.fsample;

%% Downsample & Reorder SU data
% Downsample SU data
cfg = [];
cfg.resamplefs = fs;
dataNLX = ft_resampledata(cfg, dataNLX);

% reorder
[~, b] = sort(idxElecs1);
dataNLX.label = dataNLX.label(b);
dataNLX.trial{1} = dataNLX.trial{1}(b, :);

%% Detrend and notch both datasets
% sometimes the noise is really helpful so consider commenting this part
% out
HPcutoff = 1; % high-pass filter cutoff
LPcutoff = 200; % low-pass filter cutoff
powerLineF0 = 60;
cfg = [];
cfg.bsfilter = 'yes';
cfg.bsfreq = (powerLineF0 * (1:round((LPcutoff+powerLineF0)/powerLineF0)) + [-1; 1])';
if dataNKT.fsample <= 512
    cfg.bsfreq = cfg.bsfreq(1:end-1, :); % for fs=512Hz, e.g. OS13, OS21
end
cfg.bsfiltord = 3;
cfg.hpfilter = 'yes';
cfg.hpfreq = HPcutoff;
cfg.hpfiltord = 3;
cfg.lpfilter = 'yes';
cfg.lpfreq = LPcutoff;
cfg.lpfiltord = 3;
dataNKT = ft_preprocessing(cfg, dataNKT);
dataNLX = ft_preprocessing(cfg, dataNLX);

%% Perform bipolar rereferencing
% Probably not helpful at this stage, but possible if you want it

% cfg = [];
% cfg.reref = 'yes';
% cfg.refchannel = 'all';
% cfg.refmethod = 'bipolar';
% dataNKT = ft_preprocessing(cfg, dataNKT);
% dataNLX = ft_preprocessing(cfg, dataNLX);

%% Get rid of noisy channels, if you have identidied them. 
% this bit relies on a patientInfo structure-- look at NM_infoPatients in 
% /home/knight/NM_SEEG as an example

% patientInfo = NM_infoPatients(patientCode);
% noisyChans = patientInfo.noisyChans;
%
% nChans = length(dataNKT.label);
% idxChanOk = true(nChans, 1);
% idxNoisy = false(nChans, 1);
% for idx = 1:nChans
%     labelTMP = split(dataNKT.label{idx}, '-');
%     crit1 = strcmp(regexp(labelTMP{1}, '\D+', 'once', 'match'), regexp(labelTMP{2}, '\D+', 'once', 'match')); % same elec label
%     crit2 = (str2double(regexp(labelTMP{2}, '\d+', 'once', 'match')) - str2double(regexp(labelTMP{1}, '\d+', 'match'))) == 1; % consecutive numbers
%     idxChanOk(idx) = crit1 && crit2;
%     idxNoisy(idx) = any(cellfun(@(x) ismember(x, noisyChans), labelTMP));
% end

% cfg = [];
% cfg.channel = find(idxChanOk);
% dataNKT = ft_selectdata(cfg, dataNKT);
% dataNLX = ft_selectdata(cfg, dataNLX);

%% robust scaling
% similarly unclear whether this would get rid of too much helpful noise

% NKT
N = length(dataNKT.label);
for idx = 1:N
yTMP = dataNKT.trial{1}(idx, :);
dataNKT.trial{1}(idx, :) = (yTMP - median(yTMP)) / iqr(yTMP);
end
% NLX
N = length(dataNLX.label);
for idx = 1:N
yTMP = dataNLX.trial{1}(idx, :);
dataNLX.trial{1}(idx, :) = (yTMP - median(yTMP)) / iqr(yTMP);
end

%% Change sign of dataNLX
% The +/- between the two systems tend to be reversed between the two 
% systems. You can toggle  this here in the cross correlation and viz 
% to improve fit.
dataNLX.trial{1} = -dataNLX.trial{1};

%% prep for cross correlation
L1 = length(dataNKT.time{1});
L2 = length(dataNLX.time{1});
L = max([L1 L2]);
N = length(dataNKT.label);
corrMat = zeros(N, 2*L-1);



%% Cross-correlate signal from both systems across all electrodes %%
% Pick electrode in the nlx data that has good features, is "interesting"
% in case that the electrodes do not match up between the two systems
% correlate the data not only across all possible time lags, but all
% possible electrodes. So here we pick electrode 13 in the NLK data and
% correlate it with every elec and time in the NKT data

for idx = 1:N
    [corrMat(idx, :), lags] = xcorr(dataNKT.trial{1}(idx, :), dataNLX.trial{1}(13, :));
end


%% plot lags and correlations
% look for peaks in the correlation values. I find the second plot to be
% more useful.
figure; imagesc(lags, 1:N, corrMat); colorbar;
figure; plot(lags, corrMat);




%% not sure why we do this
ratioTMP = zeros(N, 1);
for idx = 1:N
%     h(idx) = subplot(nRows, nCols, idx);
%     plot(lags, corrMat(idx, :));
%     title(sprintf('max=%.3g - std=%.3g\nratio=%.3g', max(corrMat(idx, :)), std(corrMat(idx, :)), max(corrMat(idx, :))/std(corrMat(idx, :))));
    ratioTMP(idx) = max(corrMat(idx, :))/std(corrMat(idx, :));
end

%% get the best correlation

corrMean = max(corrMat);
[~, offsetSU] = max(corrMean);

% plot it 
figure; plot(lags, corrMean); hold on; plot(lags(offsetSU), corrMean(offsetSU), 'r*');


%% redefine the trials appropriatley
if sign(lags(offsetSU)) < 0 % NLX/SU starts before NKT/clinical
    LNK = dataNKT.sampleinfo(2);
    LSU = dataNLX.sampleinfo(2)-lags(offsetSU);
    [Lnew, LnewIdx] = min([LNK LSU]);
    cfg = [];
    if LnewIdx == 1
        cfg.trl = [[1 Lnew]-lags(offsetSU) 0];
        dataNLX = ft_redefinetrial(cfg, dataNLX);
    else
        cfg.trl = [-lags(offsetSU) dataNLX.sampleinfo(2) 0];
        dataNLX = ft_redefinetrial(cfg, dataNLX);
        cfg.trl = [1 LSU 0];
        dataNKT = ft_redefinetrial(cfg, dataNKT);
    end
elseif sign(lags(offsetSU)) > 0 % NLX/SU starts after NKT/clinical
    LNK = dataNKT.sampleinfo(2)-lags(offsetSU);
    LSU = dataNLX.sampleinfo(2);
    [Lnew, LnewIdx] = min([LNK LSU]);
    cfg = [];
    if LnewIdx == 1
        cfg.trl = [1 Lnew 0];
        dataNLX = ft_redefinetrial(cfg, dataNLX);
    else
        cfg.trl = [[1 Lnew]+lags(offsetSU) 0];
        dataNKT = ft_redefinetrial(cfg, dataNKT);
    end
end


%% plot final results %%

h = zeros(2, 1);
for idxElec = 1:N
    h(1) = subplot(211);
    plot(dataNKT.trial{1}(37,:)); xlim([1 L]);
    xlim([0 11*10^5])
    h(2) = subplot(212);
    plot(dataNLX.trial{1}(4,:)); xlim([1 L]);
    linkaxes(h, 'y');
    ylim([-20 20]);
    xlim([0 11*10^5])
    input('Press ENTER to move forward...');
end


if flagFig > 0
    idx = 7;
    dataNKTcut = dataNKT.trial{1}(idx, offsetSU:offsetSU+L2-1);
    dataNLXcut = dataNLX.trial{idxTrial}(idx, :);
    dataNKTZ = (dataNKTcut-median(dataNKTcut))/(prctile(dataNKTcut, 90)-prctile(dataNKTcut, 10));
    dataNLXZ = -(dataNLXcut-median(dataNLXcut))/(prctile(dataNLXcut, 90)-prctile(dataNLXcut, 10));
    figure('Position', [218 -117 2439 727]);
    subplot(131); plot(dataNLX.time{idxTrial}, dataNKTZ); hold on; plot(dataNLX.time{idxTrial}, dataNLXZ); xlim([1 3]);
    subplot(132); plot(dataNLX.time{idxTrial}, dataNKTZ); hold on; plot(dataNLX.time{idxTrial}, dataNLXZ); xlim([150 153]);
    subplot(133); plot(dataNLX.time{idxTrial}, dataNKTZ); hold on; plot(dataNLX.time{idxTrial}, dataNLXZ); xlim([380 383]);
    figure('Position', [218 -117 2439 727]);
    subplot(311); plot(dataNLX.time{idxTrial}, dataNKTZ); title(sprintf('NK - %s - lag=%is (%i samples)', dataNKT.label{idx}, round(offsetSU/dataNKT.fsample), offsetSU));
    subplot(312); plot(dataNLX.time{idxTrial}, dataNLXZ); title(sprintf('NLX - %s', dataNLX.label{idx}));
    subplot(313); plot(dataNLX.time{idxTrial}, dataNKTZ); hold on; plot(dataNLX.time{idxTrial}, dataNLXZ);
end


%% final results
% if you are happy with the alignment and you plan to only use the NKT data
% then you just need to save the the offset and use that to redefine the
% trial moving forward so you can line up the offsets with the photodiode
% channel from the nlx data. If you want to use the nlx data, you will need
% to add to this analysis to get a more careful check to make sure each
% electrode in the nlx data finds its matching pair in the nkt data
print(cfg.trl)



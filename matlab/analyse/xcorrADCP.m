function [corr, lags, lagTime] = xcorrADCP(adcp1Data, bin1, adcp2Data, bin2, varargin)
%XCORRADCP Cross correlation between two adcp bins.
%
% Syntax:  
%       [correlation] = xcorrADCP(adcp1Data, bin1, adcp2Data, bin2)
%       [correlation] = xcorrADCP(adcp1Data, bin1, adcp2Data, bin2, Parameter, Value)
%       [correlation lags lagTime] = xcorrADCP(...)
%
% Inputs:
%
%       adcp1Data, adcp2Data  
%                       structure       OAS standard adcp data structures,
%                                       containing time series velocity data.
%                                       See help('loadADCP') for fields
%                                       reference.
%
%       bin1, bin2      [1 x 1]         Only one bin is correlated with another
%                                       at any given time. bin1 and bin2 contain
%                                       the indices of which bin to use from the
%                                       adcp1Data and adcp2Data sets.
%
% Optional Inputs:
%
%       Parameter       Value
%
%       maxlag          [1 x 2]         See help xcorr for definition of the
%                                       maxlag parameter.
%
%       option          string          See help xcorr for definition of the
%                                       option parameter.
%
%       plot            [1 x 1]         Boolean, default false
%                                       Creates a figure detailing the
%                                       correlation with time lag on the x axis
%
%       full            [1 x 1]         Boolean, default false
%                                       If true, computes full correlation
%                                       matrix of all 9 terms, otherwise just
%                                       correlation between the same velocity
%                                       components (3 column output).
%
%       trim            [1 x 1]         Boolean, default true
%                                       When true, outputs are one-sided. This
%                                       assumes that events are detected in
%                                       ADCP2 at some time AFTER they occur in
%                                       ADCP1. Outputs have the size
%                                           nC = nT + 1;
%                                       When false, the full two-sided cross
%                                       correlation is returned, so correlations
%                                       show when events happen in ADCP2 both
%                                       before and after they occur in ADCP1.
%                                       Outputs have the size
%                                           nC = 2nT + 1;
%                                       Where in both cases, nT is the lesser
%                                       of 1. the number of time points in the
%                                       data series or 2. the specified maxlag
%                                       parameter.
%
%
% Outputs:
%
%   	correlation     [nC x 3]        Correlation between the first bin and
%                                       the second, determined using xcorr(),
%                                       where nC is defined by the 'trim'
%                                       parameter above.
%                                       Columns 1 to 3 are:
%                                           u1 x u2
%                                           v1 x v2
%                                           w1 x w2
%                       or
%                       [nC x 9]        Correlation as above but including all
%                                       correlation terms. Columns 1 to 9 are:
%                                           u1 x u2
%                                           v1 x v2
%                                           w1 x w2
%                                           u1 x v2
%                                           u1 x w2
%                                           v1 x u2
%                                           v1 x w2
%                                           w1 x u2
%                                           w1 x v2
%
%       lags            [(2nT + 1) x 1] Lag indices as produced by xcorr
%       
%       lagTimes        [(2nT + 1) x 1] Lag Times from the adcp1Data time base. 
%                                       NB assumes that time based for the two
%                                       ADCPs are the same.
%
% Future Improvements:      
%
%       [1] Addition of GPUARRAY capability which is dependent on the parallel
%           processing toolbox
%
% Other m-files required:   none
% Subfunctions:             none
% Nested functions:         none
% MAT-files required:       none
%
% Author:           T. H. Clark
% Work address:     Hauser Forum
%                   3 Charles Babbage Road
%                   Cambridge
%                   CB3 0GT
% Email:            tom.clark@oceanarraysystems.com
% Website:          www.oceanarraysystems.com
%
% Created:          21 July 2014
% Revisions:        

% Set default parameters
opts.option = 'none';
opts.maxlag = [];
opts.plot   = false;
opts.full   = false;
opts.trim   = true;
opts.useGPU = false;
opts = parse_pv_pairs(opts,varargin);


%% CHECK ON INPUTS

% Check the time signals are the same


%% SIGNAL PROCESSING

% Get signal 1
u1 = adcp1Data.u(bin1,:);
v1 = adcp1Data.v(bin1,:);
w1 = adcp1Data.w(bin1,:);

% Get signal 2
u2 = adcp2Data.u(bin2,:);
v2 = adcp2Data.v(bin2,:);
w2 = adcp2Data.w(bin2,:);

% Perform correlation with the required options
if isempty(opts.maxlag)
    [corr(:,1), lags] = xcorr(u1,u2,opts.option);
    corr(:,2) = xcorr(v1,v2,opts.option);
    corr(:,3) = xcorr(w1,w2,opts.option);
    if opts.full
        corr(:,4) = xcorr(u1,v2,opts.option);
        corr(:,5) = xcorr(u1,w2,opts.option);
        corr(:,6) = xcorr(v1,u2,opts.option);
        corr(:,7) = xcorr(v1,w2,opts.option);
        corr(:,8) = xcorr(w1,u2,opts.option);
        corr(:,9) = xcorr(w1,v2,opts.option);
    end
else
    [corr(:,1), lags] = xcorr(u1,u2,opts.maxlag,opts.option);
    corr(:,2) = xcorr(v1,v2,opts.maxlag,opts.option);
    corr(:,3) = xcorr(w1,w2,opts.maxlag,opts.option);
    if opts.full
        corr(:,4) = xcorr(u1,v2,opts.maxlag,opts.option);
        corr(:,5) = xcorr(u1,w2,opts.maxlag,opts.option);
        corr(:,6) = xcorr(v1,u2,opts.maxlag,opts.option);
        corr(:,7) = xcorr(v1,w2,opts.maxlag,opts.option);
        corr(:,8) = xcorr(w1,u2,opts.maxlag,opts.option);
        corr(:,9) = xcorr(w1,v2,opts.maxlag,opts.option);
    end
end

% Get the lead/lag time
numel(adcp1Data.t)
numel(adcp2Data.t)
if numel(adcp1Data.t) < numel(adcp2Data.t)
    lagTime = [-1*fliplr(adcp1Data.t(2:lags(end))) adcp1Data.t(1:lags(end))];
    lagTime = lagTime - adcp1Data.t(1);
    disp here1
else
    lagTime = [-1*fliplr(adcp2Data.t(2:lags(end))) adcp2Data.t(1:lags(end))];
    lagTime = lagTime - adcp2Data.t(1);
    disp here2
end

size(u1)
size(u2)
size(lags)
size(corr)
size(lagTime)

% Force outputs to column vectors consistent with descriptions
lags = lags(:);
lagTime = lagTime(:);


%% TRIM CORRELATIONS

% If trim option is set, then negative time correlations are removed. That is,
% Events in ADCP2 occurring AFTER those events occurred in ADCP1 will be
% recorded - this is safe e.g. where bulk flow is moving from ADCP1 to ADCP2.
if opts.trim
    mask = lags>=0;
    corr = corr(mask,:);
    lags = lags(mask);
    lagTime = lagTime(mask);
end


%% PLOT
if opts.full && opts.plot
    
    raiseFigure(['Cross Correlation of U(ADCP 1) Bin ' num2str(bin1) ' with ADCP 2 Bin ' num2str(bin2)])
    clf
    subplot(3,1,1)
    plot(lagTime,corr(:,1))
    ylabel('u1 \cross u2')
    subplot(3,1,2)
    plot(lagTime,corr(:,4))
    ylabel('u1 \cross v2')
    subplot(3,1,3)
    plot(lagTime,corr(:,5))
    ylabel('u1 \cross w2')
    xlabel('Lag from ADCP1 Time')
    datetick('x')
    
    raiseFigure(['Cross Correlation of V(ADCP 1) Bin ' num2str(bin1) ' with ADCP 2 Bin ' num2str(bin2)])
    clf
    subplot(3,1,1)
    plot(lagTime,corr(:,6))
    ylabel('v1 \cross u2')
    subplot(3,1,2)
    plot(lagTime,corr(:,2))
    ylabel('v1 \cross v2')
    subplot(3,1,3)
    plot(lagTime,corr(:,7))
    ylabel('v1 \cross w2')
    xlabel('Lag from ADCP1 Time')
    datetick('x')
    
    raiseFigure(['Cross Correlation of W(ADCP 1) Bin ' num2str(bin1) ' with ADCP 2 Bin ' num2str(bin2)])
    clf
    subplot(3,1,1)
    plot(lagTime,corr(:,8))
    ylabel('w1 \cross u2')
    subplot(3,1,2)
    plot(lagTime,corr(:,3))
    ylabel('w1 \cross v2')
    subplot(3,1,3)
    plot(lagTime,corr(:,9))
    ylabel('w1 \cross w2')
    xlabel('Lag Time')
    datetick('x')
    
elseif opts.plot
    
    raiseFigure(['Cross Correlation of ADCP 1 Bin ' num2str(bin1) ' with ADCP 2 Bin ' num2str(bin2)])
    clf
    subplot(3,1,1)
    plot(lagTime,corr(:,1))
    ylabel('u1 \cross u2')
    subplot(3,1,2)
    plot(lagTime,corr(:,2))
    ylabel('v1 \cross v2')
    subplot(3,1,3)
    plot(lagTime,corr(:,3))
    ylabel('w1 \cross w2')
    xlabel('Lag Time')
    datetick('x')
    
end














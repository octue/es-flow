function [performancePDF, varargout] = herbie(herbPDF, evaluationFcn, varargin)
%HERBIE Monte-Carlo modeling code harness
%	Utilises probability density functions for input variables to run
%	monte-carlo simulations of arbitrary problems.
%
%	Input PDFs are described as PDFObjects and are sampled using a mersenne
%	twister approach, with seed 0
%
%   Herbie is trivially parallelisable by opening a matlabpool before execution.
%
% Syntax:
%
%                  performancePDF = herbie(a, evaluationFcn)
%                            Performs monte carlo simulation using the supplied
%                            evealuation function to evaluate performance
%
%   [performancePDF, performances] = herbie(a,evaluationFcn)
%                            Returns performances as well as the PDF
%
%   [performancePDF, performances, userResults] = herbie(a,evaluationFcn)
%                            Returns results from a customised evaluationFcn
%           
%
% Inputs:
%           a                [1 x n] PDFObject  PDFObject comprising n different
%                                               probability density functions
%                                               for n different sources of cost
%                                               (or performance, or whatever
%                                               you're modeling)
%                                               
%
%           evaluationFcn    [1 x 1] Function Handle   
%                                               Handle to a performance
%                                               evaluation function used to
%                                               build the output PDF of
%                                               performance. Example and basic
%                                               functions are given in
%                                               /monteCarlo/code_evaluationFcns
%                                               and may be added to using the 
%                                               same function template to derive
%                                               custom evaluation functions
%                                               See example below for use of an
%                                               addition function.
%
% Parameter Value Pairs (Optional Inputs):
%
%           Parameter           Value
%
%           'nEvals'            [1 x 1] double  Default 1e6
%                                               Containing integer number of
%                                               evaluations. This should be
%                                               large enough to converge the
%                                               output PDF sufficiently, small
%                                               enough that the computational
%                                               overhead is not too large
%
%           'userData'          any form        Default []
%                                               Additional data (of any form, 
%                                               usually a data structure) passed
%                                               to the evaluation function. This
%                                               is not used by herbie (hence no
%                                               constraint on form) but allows
%                                               additional data to be passed to
%                                               the evaluation function for use
%                                               when evaluating performence. For
%                                               example, one use is the
%                                               prediction of cost of energy in
%                                               a Tidal Power renewable energy
%                                               system - costs of compoonents
%                                               and operations are passed in the
%                                               input a, but the evaluation
%                                               function requires additional
%                                               data (e.g. power generation
%                                               characteristics) to perform the
%                                               performance calculation.
%
%           'showProgress'      [1 x 1] bool    Default false
%                                               Allows a progress indicator to
%                                               be displayed; useful for longer
%                                               calculations to help estimate
%                                               the finish time.
%
%           'nBins'             [1 x 1] double  Default 100
%                                               Integer number of bins in the
%                                               output PDF. Less requires less
%                                               function evaluations to
%                                               converge; more gives a finer
%                                               (better represented) output
%                                               distribution
%
% Outputs:
%
%   performancePDF              [p x 1] struct  Structure containing output
%                                               performance pdfs. p represents
%                                               the number of performance
%                                               indices returned by the
%                                               evaluation function; usually 1.
%                                               Structure has a single field:
%                             .pdf  [nBins x 2] Contains [x, PDF] values
%                                               describing the output 
%                                               probability density function
%
%   userResults             [nEvals x 1] array  An array of the additional
%                                               results returned by the
%                                               evaluation function at each
%                                               iteration.
%           
% Example:
%
%     % This example performs simple addition through a monte carlo technique
%     % Cost of 4 separate components expressed as PDF objects:
%       a(1) = PDFObject(false); % Continuous PDF
%       a(1) = setGaussianPDF(a(1),10,0.8);
%       a(2) = PDFObject(false); % Continuous PDF
%       a(2) = setGaussianPDF(a(1),15,1.4);
%       a(3) = PDFObject(true); % Discrete PDF
%       a(3) = setPDF(a(1),[2 3 4 5], [0.25 0.25 0.25 0.25]);
%       a(4) = PDFObject(false); % Continuous PDF
%       a(4) = setPiecewiseLinearPDF(a(1),[2 3 4 5], [0 1 2 0]);
%     % Run herbie using an addition evaluation function:
%       [performancePDF] = herbie(a, @herbieAddFcn)
%
% Author:                   T. H. Clark
% Work address:             Ocean Array Systems Ltd
%                           Hauser Forum
%                           3 Charles Babbage Road
%                           Cambridge
%                           CB3 0GT
% Email:                    tom.clark@oceanarraysystems.com
% Website:                  www.oceanarraysystems.com
%
% Copyright (c) 2016 Ocean Array Systems, All Rights Reserved.

% Default run parameters for herbie
herbOpts.nEvals = 1e6;
herbOpts.userData = [];
herbOpts.showProgress = false;
herbOpts.nBins = herbOpts.nEvals/1e4;

% Parse nondefault options into the structure
herbOpts    = parse_pv_pairs(herbOpts,varargin);

% Extract variables from structure to avoid additional communication overhead in
% the parfor loop
nEvals      = herbOpts.nEvals;
userData    = herbOpts.userData;
nBins       = herbOpts.nBins;

% Preallocate output variables and useful counters etc
nParams     = numel(herbPDF);
parameters  = zeros(nEvals, nParams);
userResults = cell (nEvals, 1);

% Determine input parameters for each iteration of the Monte-Carlo technique
for iParam = 1:nParams
    
    % NB samplePDF uses rand(), which in MATLAB's default implementation is
    % a mersenne twister with seed 0
    parameters(:,iParam) = samplePDF(herbPDF(iParam), nEvals);
    
end

% Preallocate performance as a single column vector. We don't know the number of
% columns a priori but when MATLAB makes the first assignment (first iteration
% of the for loop) the vector will be expanded to an array. Thus we preallocate
% the array in two steps rather than one, but we don't require a priori
% knowledge.
performance = zeros(nEvals,1);

% Run the Monte-carlo simulation nEvals times. If you have the PCT then this is
% trivially parallelised by opening a matlabpool before execution of herbie and
% changing this for to a parfor loop
for iEval = 1:herbOpts.nEvals
    
    % Current parameter set
    params = parameters(iEval,:);
    
    % Run the evaluation Function to determine the performance of this parameter
    % set
    if nargout <= 2
        [performance(iEval,:)] = evaluationFcn(params, userData);
    elseif nargout == 3
        [performance(iEval,:), userResult] = evaluationFcn(params, userData);
        % If additional results are returned by the evaluation function, they
        % must be in the form of a single element structure. This structure is
        % appended to an [nEvals x 1] userResults structure which is returned as
        % the second output argument to the command line.
        userResults(iEval) = userResult;
    end
        
end

% Number of performance variables retrieved from the evaluation function
nPerforms = size(performance,2);

% Preallocate the output structure for the required number of dimensions
performancePDF = struct('pdf', cell(nPerforms,1));

% Convert the performances to PDFs and save in the output structure
for iPerform = 1:nPerforms
    [n, xout] = hist(performance(:,iPerform),nBins);
    performancePDF(iPerform).pdf = [xout(:) n(:)./size(performance,1)];
end

% If the performances are desired, output to the command line
if nargout > 1
    varargout{1} = [parameters performance];
end

% If userResults is desired, then return it to the command line
if nargout > 2
    varargout{2} = userResults;
end


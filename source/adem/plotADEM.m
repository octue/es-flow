function [lh, fh] = plotADEM(ar, plotType, varargin)
%PLOTADEM Plot variables from an ademResults structure.
%
%       [lh, fh] = plotADEM(ar, plotType)
%       Creates a plot from data in ademResults structure ar, described by the
%       input strin plotType.
%
%       [lh, fh] = plotADEM(ar, plotType, inds)
%       Uses the specific indices into results to plot for a limited range of
%       the wall normal heights. Otherwise plots full range if possible (e.g.
%       boundary layer or Reynolds Stress profiles) or a selected range
%       (spectra)
%
%       [lh, fh] = plotADEM(..., figH)
%       Raises an existing figure and plots into it rather than creating new.
%       Does not update axes

% Parse inputs
raise = true;
inds = getInds(plotType, ar);
if nargin > 2
    if ishandle(varargin{1})
        raise = false;
        fh = varargin{1};
    else
        inds = varargin{1};
        if nargin >3
            raise = false;
            fh = varargin{2};
        end
    end
end

switch lower(plotType)
    
    case 'ux'
        
        if raise
            fh = raiseFigure('ADEM: Mean velocity profile');
            clf
        else
            fh = figure(varargin{1});
        end
        lh = plot(ar.Ux(inds), ar.eta(inds));
        
    case 'r11'
            
        if raise
            fh = raiseFigure('ADEM: Reynolds Stress R11');
            clf
        else
            fh = figure(varargin{1});
        end
        lh = plot(ar.eta(inds),ar.R(inds,1));
        
        
    case 'r12'
            
        if raise
            fh = raiseFigure('ADEM: Reynolds Stress R12');
            clf
        else
            fh = figure(varargin{1});
        end
        lh = plot(ar.eta(inds),ar.R(inds,2));

    case 'r13'
            
        if raise
            fh = raiseFigure('ADEM: Reynolds Stress R13');
            clf
        else
            fh = figure(varargin{1});
        end
        lh = plot(ar.eta(inds),ar.R(inds,3));
        
    case 'r22'
            
        if raise
            fh = raiseFigure('ADEM: Reynolds Stress R22');
            clf
        else
            fh = figure(varargin{1});
        end
        lh = plot(ar.eta(inds),ar.R(inds,4));
        
    case 'r23'
            
        if raise
            fh = raiseFigure('ADEM: Reynolds Stress R23');
            clf
        else
            fh = figure(varargin{1});
        end
        lh = plot(ar.eta(inds),ar.R(inds,5));
        
    case 'r33'
            
        if raise
            fh = raiseFigure('ADEM: Reynolds Stress R33');
            clf
        else
            fh = figure(varargin{1});
        end
        lh = plot(ar.eta(inds),ar.R(inds,6));
        
    case 'psi11'
            
        if raise
            fh = raiseFigure('ADEM: Autospectrum S11');
            clf
        else
            fh = figure(varargin{1});
        end
        lh = plot((1:size(ar.Psi,2))',ar.Psi(inds,:,1)');
        eta = ar.eta(inds);
        lam = ar.lambdaE(inds);
        for i = 1:numel(inds)
             legendCell{i} = ['\Psi (\eta = ' num2str(eta(i),'%.2f') ', \lambda_E = ' num2str(lam(i)) ')']; %#ok<AGROW>
        end
        legend(legendCell)
        
    case 'psi12'
            
        if raise
            fh = raiseFigure('ADEM: Autospectrum S12');
            clf
        else
            fh = figure(varargin{1});
        end
        
        lh = plot((1:size(ar.Psi,2))',ar.Psi(inds,:,2)');
        eta = ar.eta(inds);
        lam = ar.lambdaE(inds);
        for i = 1:numel(inds)
             legendCell{i} = ['\Psi (\eta = ' num2str(eta(i),'%.2f') ', \lambda_E = ' num2str(lam(i)) ')']; %#ok<AGROW>
        end
        legend(legendCell)
        
    case 'psi13'
            
        if raise
            fh = raiseFigure('ADEM: Autospectrum S13');
            clf
        else
            fh = figure(varargin{1});
        end
        
        lh = plot((1:size(ar.Psi,2))',ar.Psi(inds,:,3)');
        eta = ar.eta(inds);
        lam = ar.lambdaE(inds);
        for i = 1:numel(inds)
             legendCell{i} = ['\Psi (\eta = ' num2str(eta(i),'%.2f') ', \lambda_E = ' num2str(lam(i)) ')']; %#ok<AGROW>
        end
        legend(legendCell)
        
    case 'psi22'
            
        if raise
            fh = raiseFigure('ADEM: Autospectrum S22');
            clf
        else
            fh = figure(varargin{1});
        end
        
        lh = plot((1:size(ar.Psi,2))',ar.Psi(inds,:,4)');
        eta = ar.eta(inds);
        lam = ar.lambdaE(inds);
        for i = 1:numel(inds)
             legendCell{i} = ['\Psi (\eta = ' num2str(eta(i),'%.2f') ', \lambda_E = ' num2str(lam(i)) ')']; %#ok<AGROW>
        end
        legend(legendCell)
        
    case 'psi23'
            
        if raise
            fh = raiseFigure('ADEM: Autospectrum S23');
            clf
        else
            fh = figure(varargin{1});
        end
        
        lh = plot((1:size(ar.Psi,2))',ar.Psi(inds,:,5)');
        eta = ar.eta(inds);
        lam = ar.lambdaE(inds);
        for i = 1:numel(inds)
             legendCell{i} = ['\Psi (\eta = ' num2str(eta(i),'%.2f') ', \lambda_E = ' num2str(lam(i)) ')']; %#ok<AGROW>
        end
        legend(legendCell)
        
    case 'psi33'
            
        if raise
            fh = raiseFigure('ADEM: Autospectrum S33');
            clf
        else
            fh = figure(varargin{1});
        end
        
        lh = plot((1:size(ar.Psi,2))',ar.Psi(inds,:,6)');
        eta = ar.eta(inds);
        lam = ar.lambdaE(inds);
        for i = 1:numel(inds)
             legendCell{i} = ['\Psi (\eta = ' num2str(eta(i),'%.2f') ', \lambda_E = ' num2str(lam(i)) ')']; %#ok<AGROW>
        end
        legend(legendCell)
        
    otherwise
        error('Invalid Plot Type')
end
end
function inds = getInds(plotType,ar)

    switch lower(plotType(1))
        case 'p'
            % Then we plot spectra and need sensible indices!
            eta = [0.05 0.1 0.25 0.5 0.75];
            inds = interp1(ar.eta,1:numel(ar.eta),eta,'nearest',NaN);
            inds(isnan(inds)) = [];
            
        otherwise
            % We just plot all the points
            inds = 1:numel(ar.Ux);
    end
end
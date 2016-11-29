function lh = spectrumPlotvonKarman(f, uvwBar, jLi, sigma, varargin)

[asdout] = spectrumVonKarmanBasic(jLi, uvwBar, f, sigma);

% 1 x 3 boolean mask. True plots that column.
if nargin > 6
    asdout(:,~varargin{3}) = [];
end

hold on
if nargin > 5
    lh = plot(f, bsxfun(@plus,asdout,(varargin{2}(:)').^2));
else
    lh = plot(f, asdout);
end


m = max(asdout(:));

% Add frequency limits
if (nargin > 4) && ~isempty(varargin{1})
    flim = varargin{1};
    
    lh2 = plot([flim(1) flim(1)], [0 m],'r--');
    lh3 = plot([flim(2) flim(2)], [0 m],'r--');
    
    lh = [lh;lh2;lh3];
end


function lh = spectrumPlotvonKarman(f, uvwBar, jLi, varargin)

if nargin > 3
    [asdout] = spectrumVonKarmanBasic(jLi, uvwBar, f, varargin{1});
else
    [asdout] = spectrumVonKarmanBasic(jLi, uvwBar, f);
end

% 1 x 3 boolean mask. True plots that column.
if nargin > 6
    asdout(:,~varargin{4}) = [];
end

hold on
if nargin > 5
    lh = plot(f, bsxfun(@plus,asdout,(varargin{3}(:)').^2));
else
    lh = plot(f, asdout);
end


m = max(asdout(:));

% Add frequency limits
if (nargin > 4) && ~isempty(varargin{2})
    flim = varargin{2};
    
    lh2 = plot([flim(1) flim(1)], [0 m],'r--');
    lh3 = plot([flim(2) flim(2)], [0 m],'r--');
    
    lh = [lh;lh2;lh3];
end


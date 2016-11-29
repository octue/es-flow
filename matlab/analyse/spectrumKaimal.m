function [asdout] = spectrumKaimal(jLi, uvwBar, f, varargin)

% Get values out of jLi in more familiar notation:
xLu = jLi(1,1);
xLv = jLi(1,2);
xLw = jLi(1,3);

U = uvwBar(1);

if nargin>3
    sigma = varargin{1};
else
    sigma = [1 1 1];
end

f = f(:);

% KAIMAL FROM: http://www.ict-aeolus.eu/SimWindFarm/model-windfield.html
numu = 4*xLu./U;
denu = (1 + 6*f*xLu/U).^(5/3);
asdu = sigma(1)^2*numu./denu;

numv = 4*xLv./U;
denv = (1 + 6*f*xLv/U).^(5/3);
asdv = sigma(2)^2*numv./denv;

numw = 4*xLw./U;
denw = (1 + 6*f*xLw/U).^(5/3);
asdw = sigma(3)^2*numw./denw;

% Fix the NaN bug at zero frequency
asdu(isnan(asdu)) = 0;
asdv(isnan(asdv)) = 0;
asdw(isnan(asdw)) = 0;

asdout = [asdu(:) asdv(:) asdw(:)];

function [jLiOpt, sigmaOpt, sigmaNoiseOpt] = spectrumFit(method, varargin)
%SPECTRUMFIT calls spectrumFitKaimal or spectrumFitvonKarman and prints output

switch lower(method)
    case 'von karman'
        [jLiOpt, sigmaOpt, sigmaNoiseOpt] = spectrumFitvonKarman(varargin{:});
    case 'kaimal'
        [jLiOpt, sigmaOpt, sigmaNoiseOpt] = spectrumFitKaimal(varargin{:});
end

disp('________________________________________________________________________')
disp(' ')
disp(['Length and timescales, fitted using ' method ' spectrum with TRR method'])
uvwBar = varargin{3};
printScales(uvwBar, jLiOpt, sigmaOpt)

end


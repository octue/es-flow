classdef PDFObject
    %PDFObject Manages continuous and discrete probability density functions
    
    % Public Properties
    properties
        
        % Arbitrary user-modifiable description of the PDF
        Tag = '';
        
        % User data for future data exchange
        UserData
                
    end
    
    % Public properties, protected set access
    properties (SetAccess = protected)
        
        % DATA
        
        % Logical: true if discrete, false if continuous
        discreteTF
        
        % Values (column vector) at which pdf points are taken
        x
        
        % Values (column vector) of PDF at corresponding points in x
        pdf
        
    end % end protected properties
    
    
    methods
        
        % Class constructor
        function [PDFObj] = PDFObject(discrete)
            
            % Determine whether the distribution is discrete or continuous
            % and assign default uniform PDF between x = 0 and 1
            if discrete
                
                PDFObj.x   = [0; 1];
                PDFObj.pdf = [0.5; 0.5];
                PDFObj.discreteTF = true;
                
            else
                
                PDFObj.x   = [0; 1];
                PDFObj.pdf = [1; 1];
                PDFObj.discreteTF = false;
                
            end
                
        end % end function PDFObject() ( Constructor Function )
        
        % Set data
        function [PDFObj] = setPDF(PDFObj, xInput, pdfInput)
            
            % Check they're the same length
            assert(isequal(size(xInput),size(pdfInput)), 'Input x and pdf should be the same size')
            
            % Set as column vectors always
            PDFObj.x   = xInput(:);
            PDFObj.pdf = pdfInput(:);
            
            % Normalise
            sum(pdfInput)
            if PDFObj.discreteTF && (abs(sum(pdfInput)-1) > eps('single'))
                disp('PDFObject: Input discrete PDF does not sum to 1! Normalising to compensate.')
                PDFObj = normalisePDF(PDFObj);
            end
                
            
        end % end function setPDF
        
        % Set data
        function [PDFObj] = setPiecewiseLinearPDF(PDFObj, xInput, pdfInput)
            
            % Check it isn't discrete
            if PDFObj.discreteTF
                error('PDFObject:InvalidMethod','Cannot set continuous piecewise linear for a discrete-type PDF')
            end
            
            % Set as column vectors always
            PDFObj.x   = xInput(:);
            PDFObj.pdf = pdfInput(:);
            
            % Normalise
            if PDFObj.discreteTF && (sum(pdfInput)~=1)
                disp('PDFObject: Input discrete PDF does not sum to 1! Normalising to compensate.')
            end
            PDFObj = normalisePDF(PDFObj);
            
        end % end function setPiecewiseLinearPDF
        
        % Set gaussian
        function [PDFObj] = setGaussianPDF(PDFObj, mu, sigma, varargin)
            
            % Check it isn't discrete
            if PDFObj.discreteTF
                error('PDFObject:InvalidMethod','Cannot set continuous gaussian for a discrete-type PDF')
            end
            
            % If passed, varargin carries the x bounds. Default is 3 standard deviations
            if nargin > 3
                xbounds = sort(varargin{1}); % must be ascending
            else
                xbounds = [mu mu] + 3*[-sigma sigma];
            end
            
            % x points
            PDFObj.x = linspace(xbounds(1),xbounds(2),200);
        
            % Compute probability as decaying exponential away from the point
            PDFObj.pdf = exp(  ((PDFObj.x-mu).^2)./(-2*(sigma^2))  );
            
            % Normalise
            [PDFObj] = normalisePDF(PDFObj);
            
        end % end function setGaussianPDF
        
        % Set split normal
        function [PDFObj] = setSplitNormalPDF(PDFObj, mu, sigma1, sigma2, varargin)
            
            % Check it isn't discrete
            if PDFObj.discreteTF
                error('PDFObject:InvalidMethod','Cannot set split normal distribution for a discrete-type PDF')
            end
            
            % If passed, varargin carries the x bounds. Default is 3 standard deviations
            if nargin > 4
                xbounds = sort(varargin{1}); % must be ascending
            else
                xbounds = [mu mu] + 3*[-sigma1 sigma2];
            end
            
            % x points
            PDFObj.x = linspace(xbounds(1),xbounds(2),200);
        
            % Compute probability as decaying exponential away from the mean on
            % either side
            pdf1 = exp(  ((PDFObj.x-mu).^2)./(-2*(sigma1^2))  );
            pdf2 = exp(  ((PDFObj.x-mu).^2)./(-2*(sigma2^2))  );
            
            % Place the split PDF into the PDFObject
            PDFObj.pdf = pdf1;
            PDFObj.pdf(PDFObj.x>mu) = pdf2(PDFObj.x>mu);
            
            % Normalise
            [PDFObj] = normalisePDF(PDFObj);
            
        end % end function setSplitNormalPDF
        
        % Create a PDF from a piecewise linear CDF using a monte-carlo type approach
        function [PDFObj] = createFromCDF(PDFObj, xPoints, cdf, varargin)
            
            % Warning
            %warning('MATLAB:PDFObject:createFromCDF','Hacked: Function createFromCDF assumes that the CDF data input is continuous.')
            
            % Optional varargin gives number of intervals in the PDF to use. Default 200
            nSamples = 20000; % nSamples = 100*nIntervals
            if nargin>3
                nSamples = varargin{1}*100;
            end
            
            % Check that CDF(1) = 0 and CDF(end) = 1
            if (cdf(1)~=0) || (cdf(end)~=1)
                error('PDFObject:InvalidMethod','Cumulative Density Function must be bounded between 0 and 1.')
            end
            
            % Check for flats (cause singular behaviour)
            if any(cdf(1:end-1)==cdf(2:end))
                error('PDFObject:InvalidMethod','Cumulative Density Function has flat zones. Sorry, these are currently not handled.')
            end
            
            % Get the probabilities (evenly sampled)
            samples = linspace(0,1,nSamples);
            
            % Interpolate the points to give outputs
            sampledPoints = interp1(cdf,xPoints,samples,'linear');
            
            % Use hist to automatically set up bins
            [n xout] = hist(sampledPoints,nSamples/100);
            
            % Check that the object is of continuous type (required for use with
            % a continuous PDF)
            if PDFObj.discreteTF == true
            	error('PDFObject:InvalidMethod','Function createFromCDF assumes a continuous CDF. Thus, the PDFObject must be initialised as a continuous type PDF to match.')
            end
            
            % The finely discretised, discrete PDF is represented by:
            %     PDFObj.x   = xout;
            %     PDFObj.pdf = n/nSamples;
            %     PDFObj = normalisePDF(PDFObj);
            % But due to the continuous nature of the input CDF, we need a
            % continuous output PDF. Set one using a piecewise linear
            % approach (the normalisation is arbitrary based on the
            % discretisation width).
            [PDFObj] = setPiecewiseLinearPDF(PDFObj, xout, n/nSamples);
            
            
        end % end function createFromCDF    
        
        % Take N samples using the PDF
        function [xPoints] = samplePDF(PDFObj, nSamples)
            
            % Different sampling methods for discrete and continuous variables
            if PDFObj.discreteTF
                xPoints = discreteSample(PDFObj, nSamples);
            else
                xPoints = continuousSample(PDFObj, nSamples);
            end
            
        end % end function samplePDF
            
        
    end % end public methods
    
    
    
    % The following methods are hidden and can only be called by other methods
    methods (Access = protected, Hidden = true)
        
        function [PDFObj] = normalisePDF(PDFObj)
            if PDFObj.discreteTF
                PDFObj.pdf = PDFObj.pdf./sum(PDFObj.pdf);
            else
%                 PDFObj.x
%                 PDFObj.pdf
                sumPDF = trapz(PDFObj.x,PDFObj.pdf);
                disp(['    PDFObject.m:normalisePDF    Normalisation factor = ' num2str(sumPDF)])
                PDFObj.pdf = PDFObj.pdf/sumPDF;
%                 sumPDF = trapz(PDFObj.x,PDFObj.pdf)
            end
        end % end function normalisePDF
        
        function [x] = discreteSample(PDFObj, n)
            % Samples from a discrete distribution
            %       independently draws n samples (with replacement) from the 
            %       distribution specified by pdf, where pdf is a probability array 
            %       whose elements sum to 1.
            %
            %       Suppose the sample space comprises K distinct objects, then
            %       p should be an array with K elements. In the output, x(i) = k
            %       means that the k-th object is drawn at the i-th trial.

            % process p if necessary
            p = PDFObj.pdf;
            x = PDFObj.x;

            K = numel(p);
            if ~isequal(size(p), [1, K])
                p = reshape(p, [1, K]);
            end

            % construct the bins
            edges = [0, cumsum(p)];
            s = edges(end);
            if abs(s - 1) > eps
                edges = edges * (1 / s);
            end

            % draw bins
            rv = rand(1, n);
            c = histc(rv, edges);
            ce = c(end);
            c = c(1:end-1);
            c(end) = c(end) + ce;

            % extract samples
            xv = find(c);
            if numel(xv) == n  % each value is sampled at most once
                x = xv;
            else                % some values are sampled more than once
                xc = c(xv);
                d = zeros(1, n);
                dv = [xv(1), diff(xv)];
                dp = [1, 1 + cumsum(xc(1:end-1))];
                d(dp) = dv;
                x = cumsum(d);
            end

            % randomly permute the sample's order
            x = x(randperm(n));
            
            % Sample actual values not indices
            x = PDFObj.x(x);
            
        end % end function discreteSample
        
        function [x] = continuousSample(PDFObj, n)
            % Samples from a continuous distribution
            %       independently draws n samples from the 
            %       distribution specified by pdf, where pdf is a probability
            %       array whose integral across the bounds sums to 1.
            
            % Get the cumulative density function over the points
            cdf = cumtrapz(PDFObj.x,PDFObj.pdf);
            
            % Random values between 0 and 1
            probRand = rand(n,1);
            
            % Interpolate (simple linear) to convert randomly chosen
            % probabilities to x point values using the CDF
            x = interp1(cdf,PDFObj.x,probRand);
            
            
                
            
        end % end function continuousSample
        
    end
    
    
    
end % end class PDFObject






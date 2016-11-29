% Create PDFs

% Cost of 4 components into a structure
a(1) = PDFObject(false); % Continuous PDF
a(1) = setGaussianPDF(a(1),10,0.8);

a(2) = PDFObject(false); % Continuous PDF
a(2) = setGaussianPDF(a(1),15,1.4);

a(3) = PDFObject(true); % Discrete PDF
a(3) = setPDF(a(1),[2 3 4 5], [0.25 0.25 0.25 0.25]);

a(4) = PDFObject(false); % Continuous PDF
a(4) = setPiecewiseLinearPDF(a(1),[2 3 4 5], [0 1 2 0]);

% Run herbie
[performancePDF] = herbie(a, @herbieAddFcn)


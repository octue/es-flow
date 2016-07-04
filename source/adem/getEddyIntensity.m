function [J, lambda, X, Y, Z, U, V, W] = getEddyIntensity(type, varargin)
%GETEDDYINTENSITY Returns eddy intensity functions Jij, with i=1:3, j = 1:3
% These are the signatures (in terms of turbulent fluctuations) that an
% individual structure will contribute to a boundary layer.
%
% Syntax:
%       [J, lambda] = getEddyIntensity('type')
%       Returns eddy intensity functions J(i,j) for the eddy type described in
%       the input string, presently 'A', 'B1', 'B2', 'B3', 'B4' which are the
%       eddies described in Ref 1.
%
%       [J, lambda, X, Y, Z, U, V, W] = getEddyIntensity('type')
%       Outputs grids X, Y, Z, U, V, W (normalised to the eddy scale) which were
%       used to compute the eddy influence. The grids can be used to visualise
%       the influence a single eddy has on the velocity field.
%
% Inputs:
%
%       type            string      The eddy type for which to calculate J.
%                                   Current accepted values are  'A', 'B1',
%                                   'B2', 'B3', 'B4', which correspond with the
%                                   eddy types in Ref [1].
%
% Outputs:
%
%       J               [nZ x 6]    Jij contains the eddy intensity
%                                   function Jij(lambda), which is [nZ x 1] in
%                                   size. Note that the lower diagonal is not
%                                   included due to symmetry of the Reynolds
%                                   Stress Tensor. Array columns contain Jij,
%                                   ordered as:
%                                       [J11 J12 J13 J22 J23 J33];
%
%       lambda          [nZ x 1]    The remapped wall coordinate at which J is
%                                   given. Varies from 0 (at the edge of the BL)
%                                   to a point close to the wall.
%
%       X,Y,Z           [nY, nX, nZ]
%                                   Regular arrays as produced by meshgrid
%                                   containing the volume coordinate points at
%                                   which the velocity signaturee was computed.
%
%       U,V,W           [nY, nX, nZ]
%                                   Velocities in m/s induced by a single eddy
%                                   of scale z/delta = 1 and of strength gamma =
%                                   1, together with its image in the wall (z =
%                                   0), computed on the regular grid in X,Y,Z.
%
% See Also: GETREYNOLDSSTRESSES.M
%
% Future Improvements:
%
%   [1] Computation could be done using the eddy spectral function gij (eq.44
%       Perry and Marusic 1995) which is quicker but more complicated than the
%       direct integration of eqn 35 that is presently implemented.
%
%   [2] Presently the intensity functions are computed each time this is called.
%       There should be a cached lookup for particular eddy types (i.e. a .mat
%       file storing all Jij).
%
%   [3] Check why we're discretising the line filaments into small pieces -
%       probably a hangover from the old biot savart code. Eliminating (setting
%       nEl = 1) could save a lot of compute time right there (factor 51
%       improvement).
%
%   [4] Possible additional modification to take into account a free surface
%       image.
%
% References:
%
%   [1] Perry AE and Marusic I (1995) A wall-wake model for turbulent boundary
%       layers. Part 1. Extension of the attached eddy hypothesis J Fluid Mech
%       vol 298 pp 361-388
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
% Revision History:        	18 April 2015       Created from the deprecated getJ
%                                               Added the string type specifier
%                                               and reworked for new biot savart
%                                               interface.
%                           23 April 2014       Added volumetric grid and
%                                               velocity signature to
%                                               documentation
%
% Copyright (c) 2014-2015 Ocean Array Systems, All Rights Reserved.


%% If nargin = 1 we need to calculate the influence
if nargin == 1

    % EXPRESS EDDY STRUCTURE AS LINE VORTICES

    % We discretise the line vortices into smaller elements - see future
    % improvements!
    nEl = 2;

    switch lower(type)
        case 'a'

            % Define a type A eddy by its key points and reflection. See Figure 13
            % Pery and Marusic 1995.
            a = linspace(0,    1,   nEl);
            b = linspace(1,    0,   nEl);
            c = linspace(0,   -1,   nEl);
            d = linspace(-1,   0,   nEl);
            e = linspace(0,    0.8, nEl); 
            f = linspace(0.8,  0,   nEl); 
            g = linspace(0,   -0.8, nEl); 
            h = linspace(-0.8, 0,   nEl);
            L = vertcat([ a(1:end-1)'      h(1:end-1)'    a(1:end-1)'],...
                        [ b(1:end-1)'      e(1:end-1)'    b(1:end-1)'],...
                        [ a(1:end-1)'      f(1:end-1)'    c(1:end-1)'],...
                        [ b(1:end)'        g(1:end)'      d(1:end)']);

            % Start and end nodes of each vortex segment, [nNodes x 3] arrays
            startNodes  = L(1:end-1,:);
            endNodes    = L(2:end,:);

        case 'b1'

            % Define a type B1 eddy by its key points and reflection. See Figure 13
            % Pery and Marusic 1995.
            az = linspace(1,     0.5,   nEl)';
            ay = linspace(-0.15,   0,   nEl)';
            ax = linspace(0.2,     0,   nEl)';
            bz = linspace(0.5,     1,   nEl)';
            by = linspace(0,    0.15,   nEl)';
            bx = linspace(0,    0.09,   nEl)';
            arz = linspace(-1,    -0.5,   nEl)';
            ary = linspace(-0.15,    0,   nEl)';
            arx = linspace(0.2,      0,   nEl)';
            brz = linspace(-0.5,    -1,   nEl)';
            bry = linspace(0,     0.15,   nEl)';
            brx = linspace(0,     0.09,   nEl)';

            startNodes = vertcat( [ax(1:end-1)   ay(1:end-1)   az(1:end-1)],...
                                  [bx(1:end-1)   by(1:end-1)   bz(1:end-1)],...
                                  [arx(1:end-1)  ary(1:end-1)  arz(1:end-1)],...
                                  [brx(1:end-1)  bry(1:end-1)  brz(1:end-1)]);

            endNodes   = vertcat( [ax(2:end)   ay(2:end)   az(2:end)],...
                                  [bx(2:end)   by(2:end)   bz(2:end)],...
                                  [arx(2:end)  ary(2:end)  arz(2:end)],...
                                  [brx(2:end)  bry(2:end)  brz(2:end)]);

        case 'b2'

            % Define a type B2 eddy by its key points and reflection. See Figure 13
            % Pery and Marusic 1995.
            az = linspace(0.5,      1,   nEl)';
            ay = linspace(-0.15,    0,   nEl)';
            ax = linspace(0.11,   0.2,   nEl)';
            bz = linspace(1,      0.5,   nEl)';
            by = linspace(0,     0.15,   nEl)';
            bx = linspace(0.2,      0,   nEl)';
            arz = linspace(-0.5,   -1,   nEl)';
            ary = linspace(-0.15,   0,   nEl)';
            arx = linspace(0.11,  0.2,   nEl)';
            brz = linspace(-1,   -0.5,   nEl)';
            bry = linspace(0,    0.15,   nEl)';
            brx = linspace(0.2,     0,   nEl)';

            startNodes = vertcat( [ax(1:end-1)   ay(1:end-1)   az(1:end-1)],...
                                  [bx(1:end-1)   by(1:end-1)   bz(1:end-1)],...
                                  [arx(1:end-1)  ary(1:end-1)  arz(1:end-1)],...
                                  [brx(1:end-1)  bry(1:end-1)  brz(1:end-1)]);

            endNodes   = vertcat( [ax(2:end)   ay(2:end)   az(2:end)],...
                                  [bx(2:end)   by(2:end)   bz(2:end)],...
                                  [arx(2:end)  ary(2:end)  arz(2:end)],...
                                  [brx(2:end)  bry(2:end)  brz(2:end)]);

        case 'b3'

            % Define a type B3 eddy by its key points and reflection. See Figure 13
            % Pery and Marusic 1995.
            az = linspace(1,       0.5,   nEl)';
            ay = linspace(-0.15,     0,   nEl)';
            ax = linspace(0.09,      0,   nEl)';
            bz = linspace(0.5,       1,   nEl)';
            by = linspace(0,      0.15,   nEl)';
            bx = linspace(0,       0.2,   nEl)';
            arz = linspace(-1,     -0.5,   nEl)';
            ary = linspace(-0.15,     0,   nEl)';
            arx = linspace(0.09,      0,   nEl)';
            brz = linspace(-0.5,     -1,   nEl)';
            bry = linspace(0,      0.15,   nEl)';
            brx = linspace(0,       0.2,   nEl)';

            startNodes = vertcat( [ax(1:end-1)   ay(1:end-1)   az(1:end-1)],...
                                  [bx(1:end-1)   by(1:end-1)   bz(1:end-1)],...
                                  [arx(1:end-1)  ary(1:end-1)  arz(1:end-1)],...
                                  [brx(1:end-1)  bry(1:end-1)  brz(1:end-1)]);

            endNodes   = vertcat( [ax(2:end)   ay(2:end)   az(2:end)],...
                                  [bx(2:end)   by(2:end)   bz(2:end)],...
                                  [arx(2:end)  ary(2:end)  arz(2:end)],...
                                  [brx(2:end)  bry(2:end)  brz(2:end)]);

        case 'b4'

            % Define a type B4 eddy by its key points and reflection. See Figure 13
            % Pery and Marusic 1995.
            az = linspace(0.5,      1,   nEl)';
            ay = linspace(-0.15,    0,   nEl)';
            ax = linspace(0,      0.2,   nEl)';
            bz = linspace(1,      0.5,   nEl)';
            by = linspace(0,     0.15,   nEl)';
            bx = linspace(0.2,   0.11,   nEl)';
            arz = linspace(-0.5,    -1,   nEl)';
            ary = linspace(-0.15,    0,   nEl)';
            arx = linspace(0,      0.2,   nEl)';
            brz = linspace(-1,    -0.5,   nEl)';
            bry = linspace(0,     0.15,   nEl)';
            brx = linspace(0.2,   0.11,   nEl)';
            startNodes = vertcat( [ax(1:end-1)   ay(1:end-1)   az(1:end-1)],...
                                  [bx(1:end-1)   by(1:end-1)   bz(1:end-1)],...
                                  [arx(1:end-1)  ary(1:end-1)  arz(1:end-1)],...
                                  [brx(1:end-1)  bry(1:end-1)  brz(1:end-1)]);

            endNodes   = vertcat( [ax(2:end)   ay(2:end)   az(2:end)],...
                                  [bx(2:end)   by(2:end)   bz(2:end)],...
                                  [arx(2:end)  ary(2:end)  arz(2:end)],...
                                  [brx(2:end)  bry(2:end)  brz(2:end)]);

        otherwise

            error('MATLAB:InvalidInput','Unrecognised eddy type string input. See help(''getEddyIntensity'') for valid eddy types')

    end



    % PLOT EDDY STRUCTURE SHAPE

    raiseFigure(['Type ' type ' Eddy Structure'])
    clf
    subplot(1,3,1)
    % plot3(startNodes(:,1), startNodes(:,2), startNodes(:,3))
    plot3([startNodes(:,1)'; endNodes(:,1)'],[startNodes(:,2)'; endNodes(:,2)'],[startNodes(:,3)'; endNodes(:,3)'],'k-')
    axis equal
    xlabel('x/\delta')
    ylabel('y/\delta')
    zlabel('z/\delta')



    % SETUP INFLUENCE DOMAIN

    dx = 0.010;
    dy = 0.010;
    dz = 0.020;
    xVec = 0:dx:4;
    xVec = [fliplr(-1*xVec(1,2:end)) xVec];
    yVec = 0:dy:2;
    yVec = [fliplr(-1*yVec(1,2:end)) yVec];
    zVec = 0:dz:1.5;

    % Size of domain
    nX = numel(xVec);
    nY = numel(yVec);
    nZ = numel(zVec);

    U = zeros(nY, nX, nZ);
    V = zeros(nY, nX, nZ);
    W = zeros(nY, nX, nZ);

    % Volume Mesh. Note we produce with meshgrid rather than ndgrid for
    % compatibility with MATLAB's volume plotting facilities - helps us to visualise
    % the eddy signature.
    [X,Y,Z] = meshgrid(xVec,yVec,zVec);



    % DETERMINE VELOCITY INFLUENCE OF THE STRUCTURE

    % Squares of the core radii, [nNodes x 1] array
    rcEffSqd    = (0.05^2)*ones(size(startNodes,1),1);

    % Vortex strengths, [nNodes x 1] array
    gamma       = ones(size(startNodes,1),1);

    % Grid locations on which to compute influence
    locations   = [X(:), Y(:), Z(:)];

    % Type of core model to use (prevents singularities)
    CoreModel   = 2;

    % Type of computation to use (currently only 2, which uses the GPU, is allowed)
    RunMode = 2;

    % Cutoff radius (10 deltas - encompasses the whole domain for these purposes)
    CutOff      = 10;

    % GPU optimisation parameters for if a GPU is present
    GPUParams   = [1024 1 0];

    % Do the Biot Savart calculationdisp(['bsTest1: Using file ' which('biotSavart')])
    disp( 'getEddyIntensity: Running biotSavart code with parameters:' )
    disp(['       Eddy Type: ' type                     ]) 
    disp(['       CoreModel: ' num2str(CoreModel)       ]) 
    disp(['          CutOff: ' num2str(CutOff)          ])
    disp(['         RunMode: ' num2str(RunMode)         ])
    disp(['       GPUParams: [' num2str(GPUParams) ']'  ])
    disp(['        MEX File: ' which('biotSavart')      ])
    [uind, ~] = biotSavart(startNodes', endNodes', locations', gamma, rcEffSqd, CoreModel, CutOff, RunMode, GPUParams);

    % If we wish to save a test case for the biot savart code (it's a useful one
    % becuse of the fine and regular grid - NaNs pop up where the line vortex
    % elements are collinear with grid points. These NaNs are checked for and
    % handled in the BS code without upsetting the results, but it'd be good if they
    % didn't crop up in the first place...
    % save('bsTest4.mat','locations','startNodes','endNodes','gamma','rcEffSqd','CoreModel','CutOff','RunMode','GPUParams')

    U(:) = uind(1,:);
    V(:) = uind(2,:);
    W(:) = uind(3,:);

else
    [X, Y, Z, U, V, W] = deal(varargin{:});
    %     xVec = X(1,:,1); xVec = xVec(:);
    %     yVec = Y(:,1,1); yVec = yVec(:);
    zVec = Z(1,1,:); zVec = zVec(:);
end


%% COMPUTE AND PLOT REYNOLDS STRESSES

V1V1 = U.*U;
V1V2 = U.*V;
V1V3 = U.*W;
V2V2 = V.*V;
V2V3 = V.*W;
V3V3 = W.*W;

U0 = 1/(2*pi);

% Compute I
I11 = trapz(Y(:,1,1), trapz(X(1,:,1),V1V1./U0^2,2),1);
I12 = trapz(Y(:,1,1), trapz(X(1,:,1),V1V2./U0^2,2),1);
I13 = trapz(Y(:,1,1), trapz(X(1,:,1),V1V3./U0^2,2),1);
I22 = trapz(Y(:,1,1), trapz(X(1,:,1),V2V2./U0^2,2),1);
I23 = trapz(Y(:,1,1), trapz(X(1,:,1),V2V3./U0^2,2),1);
I33 = trapz(Y(:,1,1), trapz(X(1,:,1),V3V3./U0^2,2),1);

% And infer J which is equal; just remapped to a different space
% Compute the z/delta spacing from lambda the logarithmic spacing
% lambda = linspace(-10, log(1/0.0001), 1000)';
lambda = linspace(0,log(1000),50);
newZVec = 1./exp(lambda);
J11 = interp1(zVec(:),I11(:),newZVec(:),'pchip',NaN);
J12 = interp1(zVec(:),I12(:),newZVec(:),'pchip',NaN);
J13 = interp1(zVec(:),I13(:),newZVec(:),'pchip',NaN);
J22 = interp1(zVec(:),I22(:),newZVec(:),'pchip',NaN);
J23 = interp1(zVec(:),I23(:),newZVec(:),'pchip',NaN);
J33 = interp1(zVec(:),I33(:),newZVec(:),'pchip',NaN);

% Reorganise output to a structure
J = [J11(:) J12(:) J13(:) J22(:) J23(:) J33(:)];

% Eliminate values that were out of the bounds of our computed domain
mask = any(isnan(J),2);
J = J(~mask,:);
lambda = lambda(~mask);



%% PLOT EDDY INTENSITY FUNCTIONS

raiseFigure(['Type ' type ' Eddy Structure'])

subplot(1,3,2)
plot(zVec(:),[I11(:) I12(:) I13(:) I22(:) I23(:) I33(:)])
xlabel('z/\delta')
ylabel('Eddy Intensity Function I')
legend({'I_1_1';'I_1_2';'I_1_3';'I_2_2';'I_3_3'})

subplot(1,3,3)
plot(lambda,J)
xlabel('\lambda')
ylabel('Eddy Intensity Function J')
legend({'J_1_1';'J_1_2';'J_1_3';'J_2_2';'J_2_3';'J_3_3'})








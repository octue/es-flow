function [tsVec, streams, upstreams, distances, velocities, depths] = streamEulerian(s, seeds, receivers, monitors, t)
%STREAMEULERIAN Finds distances between instantaneous streamlines and points.
% Optionally interpolates velocity and depth fields to find u,v, depth at input monitor points.

      

% Ascertain streamlines thru seed locations
%   Because timesteps are taken 10 minutes apart, a lagrangian approach (where
%   particle locations are seeded at one instant then tracked over time) is
%   inappropriate, because the timestep beween the fields is too large to be of
%   use in the integration.
%   Thus, we adopt an eulerian approach. This assumes that a particle dropped in
%   the flow at a seed point will follow the instantaneous streamline passing
%   through the seed point at that moment in time.



% Counters
nSources   = size(seeds,1);
nReceivers = size(receivers,1);
nTimeSteps = size(s.u,3);
nMonitors  = size(monitors,1);

% Do one or many timesteps
if nargin > 4
    tsVec = interp1(s.t,1:nTimeSteps,t,'nearest');
    nTimeSteps = 1;
else
    % Do for all timesteps
    tsVec = 1:nTimeSteps;
end

% Preallocate output arrays
distances = zeros(nSources, nReceivers, nTimeSteps);
streams = cell(nSources,nTimeSteps);
upstreams = cell(nSources,nTimeSteps);
velocities = zeros(nMonitors, 2, nTimeSteps);
depths = zeros(nMonitors, nTimeSteps);

% For each timestep
ctr = 0;
for timeStep = tsVec;
    ctr = ctr+1;
    
    % Compute the Eulerian streamlines
    streamsTemp   = stream2(s.x, s.y, s.u(:,:,timeStep), s.v(:,:,timeStep), seeds(:,1), seeds(:,2));
    upstreamsTemp = stream2(s.x, s.y, -1*s.u(:,:,timeStep), -1*s.v(:,:,timeStep), seeds(:,1), seeds(:,2));
      
    % Save the streamline data (streamlines are in cells)
    streams(:,ctr)   = streamsTemp(:);
    upstreams(:,ctr) = upstreamsTemp(:);
    
    if ~isempty(receivers)
        % Compute minimum distance between streamline sources and receiver locations
        dist = zeros(nSources, nReceivers);
        for sCtr = 1:nSources
            for rCtr = 1:nReceivers
                dist(sCtr,rCtr) = minDistance(streamsTemp{sCtr}, receivers(rCtr,:));
            end
        end
    % Assign the matrix for storage at this timestep
    distances(:,:,ctr) = dist;
    
    end
    if ~isempty(monitors)
        % Lookup velocities at monitor points for the current timestep
        velocities(:,1,ctr) = interp2(s.x,s.y,s.u(:,:,timeStep),monitors(:,1),monitors(:,2),'spline');
        velocities(:,2,ctr) = interp2(s.x,s.y,s.v(:,:,timeStep),monitors(:,1),monitors(:,2),'spline');

        % Lookup depths at monitor points for the current timestep
        depths(:,ctr) = interp2(s.x,s.y,s.depth(:,:,timeStep), monitors(:,1), monitors(:,2),'spline');
    end
end

end % end Main Function


function dist = minDistance(XY, point)
%MINDISTANCE Computes the minimum distance between a line and a point. The line
% can be defined by one or more many straight segments. Vector projection onto
% infinitely long lines collinear with each segment is used. Bounds are applied
% where the projection of the point onto the line is outside the bounds of the
% line segment(s).

% Line equation x = a + lamda*n where a is a point on the line, n is the
% normalised unit vector of the line.

% Unit vectors for each line segment
magLine = sqrt(sum(diff(XY).^2,2));
normLine = diff(XY)./repmat(magLine, [1 2]);

% Vector from point to line segment start location
aMinusp = bsxfun(@minus,XY(1:end-1,:),point);

% Projection of (a-p) onto the line
projection = repmat(dot(aMinusp,normLine,2), [1 2]).*normLine;

% If Projection is in the same direction as normline, then the point is outside
% the bounds and a is the closest point. If projection is in the opposite
% direction to normline but larger thanthe magnitude of the line segment, then
% it is out of bounds and the endpoint is the closest. Otherwise the projection
% of the point onto the line is within the bounds of the segment, and the
% orthogonal vector will give the minimum distance

% Orthogonal vector from point to line
orthVector = aMinusp - projection;

% Norm of the orthogonal vector gives closest point between infinite line and
% the point
dist = sqrt(sum(orthVector.^2,2));

% Out of bounds before the line start point a
oob_a = sign(projection(:,1))==sign(normLine(:,1));
adist = sqrt(sum(aMinusp.^2,2));
dist(oob_a) = adist(oob_a);

% Out of bounds beyond the line end point b
oob_b = (~oob_a) & (sqrt(sum(projection.^2,2)) > magLine);
bMinusp = bsxfun(@minus,XY(2:end,:),point);
bdist = sqrt(sum(bMinusp.^2,2));
dist(oob_b) = bdist(oob_b);

% There can be many line segments; we're looking for the closest point on any of
% them
dist = min(dist);

end
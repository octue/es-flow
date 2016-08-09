function [unitPrice] = herbieMECFcn(params, userData)
%HERBIEMECFCN Performs Tidal Energy Cost Modeling

% Note that this function DOES NOT attempt to resolve maintenance
%                               dependencies, nor does it optimise maintenance
%                               scheduling. You are left to sort out maintenance
%                               programming for yourself!

% Get structures out of the user data array
fixCapex = userData.fixCapex;
varCapex = userData.varCapex;
fixOpex  = userData.fixOpex;
varOpex  = userData.varOpex;

% Establish numbers of items of each type
nFixCapex = size(fixCapex.itemNames,1);
nVarCapex = size(varCapex.itemNames,1);
nFixOpex  = size(fixOpex.itemNames, 1);
nVarOpex  = size(varOpex.itemNames, 1);

% Split the parameters array into the contituent cost structures
fixCapex.cost = params(1:nFixCapex);
params(1:nFixCapex) = [];
varCapex.cost = params(1:nVarCapex);
params(1:nVarCapex) = [];
fixOpex.cost = params(1:nFixOpex);
params(1:nFixOpex) = [];
varOpex.cost = params(1:nVarOpex);

% Get ppkwh inputs out of userData
unitPrice = ppkwh(  userData.discountRate, ...
                    userData.capacity, ...
                    userData.production, ...
                    fixCapex, varCapex, fixOpex, varOpex);

end


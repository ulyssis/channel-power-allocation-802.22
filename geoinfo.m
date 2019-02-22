% generate location of DTT, WBS and terminals, and the attenuation among
% them
function [posSU, posET, posTVContor, Gtilde, GtildeETsSUs, GtildeAll] = geoinfo(n, m, nET, lengthSide, coverage, SUcellRadius, pathlossfactor, s)


    [posSU, posET]= SUPositionInGrid(n, nET, lengthSide, coverage);        % Randomly locate the SUs&ETs within a lengthSide x lengthSide area
%                 [posSU]= SUPositionInGrid2(n, lengthSide);        % Randomly locate the SUs within a 10 x 10 area
    [posTVContor]= PUContorPosition(m, lengthSide);     % Locate the contors of PUs
    [posTV]= PUPositionFixed(m, lengthSide);     % Locate Primary TV stations which cause interference to SUs

    %%% the pathloss matrix can be integreted into a big one,
    %%% but will not be flexible to cancel shadowing in a
    %%% particular part.
    % pathloss(atinuation+shadowing) among SUs
    d = dist(posSU);    
    G = ones(n,n)./ (d - SUcellRadius + diag(ones(n, 1))).^pathlossfactor;
    shadow = normrnd (0, s, n, n);
    G = G.*10.^(shadow/10) ;
    Gtilde = G - diag(diag(G));     % path loss between secondary base station to the nearest edge of other base stations

    % pathloss(atinuation+shadowing) among ETs and SUs
    posETsSUs = [posET, posSU];
    d = dist(posETsSUs);    % get distance matrix for SUs
    GETs = ones(n+n*nET, n+n*nET)./ (d  + diag(ones(n+n*nET, 1))).^pathlossfactor;
    shadow = normrnd (0, s, n+n*nET, n+n*nET);
    GETs = GETs.*10.^(shadow/10);
    GtildeETsSUs = GETs - diag(diag(GETs));     % path loss between any ET of one WBS and any other WBS     


    % pathloss(atinuation+shadowing) among SUs and TV contours
    d = dist([posSU posTVContor posTV]);      % get channel gain matrix for SUs + TVContors + TVs
    GAll = ones(n+m+m, n+m+m)./ (d - SUcellRadius + diag(ones(n+m+m, 1))).^pathlossfactor;
%                 shadow = normrnd (0, s, n+m+m, n+m+m);
%                 GAll = GAll.*10.^(shadow/10);
    GtildeAll = GAll - diag(diag(GAll));





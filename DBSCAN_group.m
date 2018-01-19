function [clusterID,pntType,numCluster,clusterSize] = DBSCAN_group(pntNN,pntScore,critScore,varargin)
% DBSCAN_group implements point clustering based on the principle of 
% DBSCAN (Density-Based Spatial Clustering of Applications with Noise). 
% DBSCAN is a data clustering algorithm proposed by Martin Ester, 
% Hans-Peter Kriegel, Jörg Sander and Xiaowei Xu in 1996.
%
%   INPUTS:
%   pntNN: cell array of vectors; contains for each point the list indices of all potential links  
%   pntScore: vector; density score for each point
%   critScore: threshold for core point class
%
%   written by
%   C.P.Richter
%   Division of Biophysics / Group J.Piehler
%   University of Osnabrueck
%
%   modified 18.05.2015: fixed a bug that allowed the clusters to propagate through border points

ip = inputParser;
ip.KeepUnmatched = true;
addRequired(ip,'pntNN')
addRequired(ip,'pntScore')
addRequired(ip,'critScore',@(x)isscalar(x) && x > 0)
addParamValue(ip,'pntOrder', 'ascend', @(x)ischar(x))
addParamValue(ip,'verbose', false, @(x)islogical(x))
parse(ip,pntNN,pntScore,critScore,varargin{:});

pntOrder = ip.Results.pntOrder;   
verbose = ip.Results.verbose;

%%
numPnts = numel(pntScore);
switch pntOrder
    case 'ascend'
        pntList = 1:numPnts;
    case 'descend'
        pntList = numPnts:-1:1;
end %switch

pntType = (pntScore >= critScore); %is core point

clusterID = nan(numPnts,1); %unclassified
clusterIdx = 0;
for pntIdx = pntList
    if pntType(pntIdx) == 1 && isnan(clusterID(pntIdx))
        %initialize new cluster
        clusterIdx = clusterIdx+1;
        
        %point is put into respective cluster
        clusterID(pntIdx) = clusterIdx;
        
        isConn = pntNN{pntIdx};
        take = isnan(clusterID(isConn)); %only use those points not already classified
        while any(take)
            clusterID(isConn(take)) = clusterIdx;
            %find all connected points (density-reacheability)    
            isConn_ = unique(horzcat(pntNN{isConn(take)}));
            isConn = isConn_(not(ismembc(isConn_,isConn(take))));
            %                         isConn = setdiff(horzcat(pntNN{isConn(take)}),isConn(take));
            take = isnan(clusterID(isConn)) & (pntType(isConn) == 1); %make sure clusters propagate only through core points
        end %while
    end %if
end %for
clusterID = clusterID + 1;
isNoise = isnan(clusterID);
clusterID(isNoise) = 1; %assign noise points to the same cluster

pntType = double(pntType);
pntType(isNoise) = -1; %is noise point

numCluster = clusterIdx; %(= #, noise cluster excluded)
clusterSize = accumarray(clusterID,1);

%%
if verbose
    %%
    [f,xbin] = hist_fd(clusterSize(2:end));
    figure; hold on
    plot(xbin,f,'k.')
    xlabel('Clustersize')
    axis tight
    box on
    
    %%
    [f,xbin] = hist(pntType,[-1 0 1]);
    figure; hold on
    bar(xbin,f,'hist')
    xlabel('Pointtype')
    axis tight
    box on
end %if
end %fun
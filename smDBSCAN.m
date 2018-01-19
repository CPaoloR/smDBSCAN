function [clusterID,pntType] = smDBSCAN(SML,searchRad,varargin)
% smDBSCAN applies the DBSCAN principle of density-reachability to group 
% subsequent observations of the same immobile single-molecule over time.
%
%   INPUTS:
%   SML: structure; stores the single-molecule localization information
%   SML.t: vector; timepoint (image frame) at which the molecule has been localized
%   SML.i: vector; i-position of the molecule
%   SML.j: vector; j-position of the molecule
%   searchRad: scalar; defines the max. distance between observations to be considered as potential links
%
%   written by
%   C.P.Richter
%   Division of Biophysics / Group J.Piehler
%   University of Osnabrueck
%
%% 
ip = inputParser;
ip.KeepUnmatched = true;
addRequired(ip,'SML')
addRequired(ip,'searchRad', @(x)isscalar(x) & x > 0)
addParamValue(ip,'xFactor',0.75) % times 100% of the current time-window
addParamValue(ip,'minT', 30) %[frames]
addParamValue(ip,'maxT', inf) %[frames]
addParamValue(ip,'numT', 15)
addParamValue(ip,'scaleT', 'log10') % or linear
addParamValue(ip,'verbose', false, @(x)islogical(x))
parse(ip,SML,searchRad,varargin{:});

xFactor = ip.Results.xFactor;
minT = ip.Results.minT;
maxT = ip.Results.maxT;
numT = ip.Results.numT;
scaleT = ip.Results.scaleT;
verbose = ip.Results.verbose;

if isinf(maxT)
    maxT = range(SML.t);
end %if
switch scaleT
    case {'log10','lg'}
        T = unique(round(logspace(log10(minT),log10(maxT),numT))); %[frames] temporal search windows used in the iterative filtering process
    case {'linear','lin'}
        T = unique(round(linspace(minT,maxT,numT)));
end %switch

%% CLUSTERING
%% iterative segmentation of the clusters
take = 1:numel(SML.t); %we start with all points in the game
clusterID = ones(numel(SML.t),1); %initialize as noise
pntType = -1*ones(numel(SML.t),1); %initialize as noise

%starting from the long lasting clusters as these will be picked up with higher fidelity
%(at smaller T, the chance for false identification is increased due to the less stringent requirements)
for idxT = numT:-1:1
    numCluster = max(clusterID) - 1; % ID = 1 being the noise cluster
    
    %search for pot. links (observations that potentially stem from the same emitter)
    pntNN = DBSCAN_pot_link([SML.i(take) SML.j(take)],searchRad(1));
    
    critScore = xFactor*T(idxT);
    %apply forward and backward clustering with DBSCAN
    [clusterID(take),pntType(take)] = ...
        DBSCAN_fwd_rev_cluster(...
        get_cross_field_value(SML,take),...
        pntNN,T(idxT),...
        critScore,...
        'idOffset',numCluster);
    
    %% start next iteration on the set of unclassified points
    take = find(pntType == -1);
end %for

clear pntNN

%%
if verbose
    show_me_my_cluster(SML,clusterID,pntType)
end %if
end %fun
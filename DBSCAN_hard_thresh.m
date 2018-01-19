function pntNN = DBSCAN_hard_thresh(SML,pntNN,timeWin,varargin)
% DBSCAN_hard_thresh reduces the list of potential links to those 
% that reside within the expected time window around each point
%
%   INPUTS:
%   SML: structure; stores the single-molecule localization information
%   SML.t: vector; timepoint (image frame) at which the molecule has been localized
%   SML.i: vector; i-position of the molecule
%   SML.j: vector; j-position of the molecule
%   pntNN: cell array of vectors; contains for each point the list indices of all potential links  
%   timeWin: vector; defines the time window around each point
%
%   written by
%   C.P.Richter
%   Division of Biophysics / Group J.Piehler
%   University of Osnabrueck
%
%   modified 19.05.2015: fixed a bug potentially leading to a false time-coordinate for the query point

ip = inputParser;
ip.KeepUnmatched = true;
addRequired(ip,'SML',@iscolumn)
addRequired(ip,'pntNN',@iscell)
addRequired(ip,'timeWin')

parse(ip,SML,pntNN,timeWin,varargin{:});

%%
for pntIdx = numel(pntNN):-1:1
    take = ismembc(SML.t(pntNN{pntIdx}),SML.t(pntIdx)+timeWin);
    
    pntNN{pntIdx} = pntNN{pntIdx}(take);
    
    if numel(pntNN{pntIdx}) > 1
        %check if there are multiple observations within 1 frame
        numObsT = accumarray(SML.t(pntNN{pntIdx}),1);
        obsT = find(numObsT);        
        simultanObsT = (nonzeros(numObsT) > 1);
        if any(simultanObsT)
            %discard all multiple observation (only the NN remains)
            for idxSimultanObsT = reshape(find(simultanObsT),1,[])
                isSimultanObsT = find(SML.t(pntNN{pntIdx}) == obsT(idxSimultanObsT));
                
                pntDist = sqrt((SML.i(pntNN{pntIdx}(isSimultanObsT))-SML.i(pntIdx)).^2 + ...
                    (SML.j(pntNN{pntIdx}(isSimultanObsT))-SML.j(pntIdx)).^2);                
                [~,idxNN] = min(pntDist);
                isSimultanObsT(idxNN) = []; %remove the nearest-neighbor from the list
                pntNN{pntIdx}(isSimultanObsT) = []; %remove all other oberservations from the neighborhood list
            end %for
        end %if
    end %if
end %for
end %fun
function [pntNN,pntDist] = DBSCAN_pot_link(X,searchRad,varargin)
% DBSCAN_POT_LINK Find points that potentially belong into one group
%   INPUTS:
%   X -> N x D matrix; encoding the position of each observation in space
% searchRad -> scalar; defines the max. distance between observations to be considered as potential links

%written by
%C.P.Richter
%Division of Biophysics / Group J.Piehler
%University of Osnabrueck

%modified 08.05.2015: added conversion to uint32 to save 50%

%%
ip = inputParser;
ip.KeepUnmatched = true;
addRequired(ip,'X')
addRequired(ip,'searchRad')
addParamValue(ip,'verbose', false, @(x)islogical(x))
parse(ip,X,searchRad,varargin{:});

verbose = ip.Results.verbose;

%% calculate the nearest-neighbor relationship
if all(searchRad == searchRad(1))
    objKdTree = KDTreeSearcher(X);
    if nargout == 1
        pntNN = rangesearch(objKdTree,X,searchRad(1));
    else
        [pntNN,pntDist] = rangesearch(objKdTree,X,searchRad(1));
    end %if
else %normalize each dimension so a distance 1 is equal that search radius (ball/ellipsoid search)
    X = bsxfun(@times,X,1./searchRad);
    objKdTree = KDTreeSearcher(X);
    if nargout == 1
        pntNN = rangesearch(objKdTree,X,1);
    else
        [pntNN,pntDist] = rangesearch(objKdTree,X,1);
    end %if
end

pntNN = cellfun(@uint32,pntNN,'un',0); %saves 50% RAM compared to double precision

%%
if verbose
    %%
    numPotLink = cellfun('size',pntNN,2);
    if not(isempty(numPotLink))
        [f,x] = ecdf(numPotLink);
        
        hFig = figure('Color','w'); hold on
        plot(x,f,'linewidth',2,'color','k')
        xlabel('# pot. links','FontSize',20)
        ylabel('CDF','FontSize',20)
        axis tight
        box on
        set(gca(hFig),'FontSize',20)
    end %if
end %if
end %fun
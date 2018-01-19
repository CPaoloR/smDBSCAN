function [hAx,hFig] = show_me_my_cluster(SML,clusterID,pntType,varargin)
% show_me_my_cluster visualizes the results of the clustering process
%
%   INPUTS:
%   SML: structure; stores the single-molecule localization information
%   SML.t: vector; timepoint (image frame) at which the molecule has been localized
%   SML.i: vector; i-position of the molecule
%   SML.j: vector; j-position of the molecule
%
%   written by
%   C.P.Richter
%   Division of Biophysics / Group J.Piehler
%   University of Osnabrueck

%%
ip = inputParser;
ip.KeepUnmatched = true;
addParamValue(ip,'Title', '')
parse(ip,varargin{:});

%% visualization
hFig = figure('Color','w');
hold on

plot3(SML.j,SML.i,SML.t,'k.')
%noise
take = (pntType == -1) | (pntType == 0);
plot3(SML.j(take),SML.i(take),SML.t(take),'k.')
clusterIDs = unique(clusterID);
numCluster = numel(clusterIDs)-1;

if numCluster > 0
    clusterColor = hsv2rgb([colvec(randperm(numCluster))/...
        numCluster,ones(numCluster,2)]);
    for idxCluster = 2:numCluster
        %     if idxCluster == 1
        %         %skip noise cluster
        %         continue
        %     end
        %     if clusterSize(idxCluster) > 10
        %     clusterColor = hsv2rgb([rand(1),1,1]);
        
        %     take = (clusterID == idxCluster) & (pntType == 1);
        take = (clusterID == clusterIDs(idxCluster));
        if any(take)
            plot3(SML.j(take),SML.i(take),SML.t(take),'o',...
                'color',clusterColor(idxCluster,:),'markersize',6)
            
            %%mark caps
            idxMinCap = min(SML.t(take));
            plot3(mean(SML.j(take)),mean(SML.i(take)),idxMinCap,...
                'square','color',clusterColor(idxCluster,:),'markersize',20,'linewidth',2)
            idxMaxCap = max(SML.t(take));
            plot3(mean(SML.j(take)),mean(SML.i(take)),idxMaxCap,...
                'square','color',clusterColor(idxCluster,:),'markersize',20,'linewidth',2)
            
            %%mark border points
            take = (clusterID == idxCluster) & (pntType == 0);
            plot3(SML.j(take),SML.i(take),SML.t(take),'o',...
                'markersize',6,'color',clusterColor(idxCluster,:))
        else
            fprintf('Problem with ClusterID = %d\n',idxCluster)
        end %if
    end %for
end %if
xScale = [floor(min(SML.j)),ceil(max(SML.j))];
yScale = [floor(min(SML.i)),ceil(max(SML.i))];
zScale = [floor(min(SML.t)),ceil(max(SML.t))];

hAx = gca(hFig);
title(ip.Results.Title)

% ylabel('y [px]','Fontsize',20); ...
%     xlabel('x [px]','Fontsize',20); ...
box on
set(hAx,'units','normalized','outerposition',[0.1 0.1 0.8 0.8],...
    'dataaspectratio',[1 1 range(SML.t)./max(range(SML.i),range(SML.j))],'FontSize',20,...
    'Xlim',xScale,'Xtick',[xScale(1) mean(xScale) xScale(2)],...
    'Xticklabel',{num2str(xScale(1),'%d'),'x [px]',num2str(xScale(2),'%d')},...
    'Ylim',yScale,'Ytick',[yScale(1) mean(yScale) yScale(2)],...
    'Yticklabel',{num2str(yScale(1),'%d'),'y [px]',num2str(yScale(2),'%d')},...
    'Zlim',zScale)
axis vis3d ij
zlabel('Time [frame]','Fontsize',20)
end
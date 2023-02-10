function [hb_pruned] = prune_channels_nan(hb,ScansQuality)
% This functions prunes channels based on QT-nirs results. 
% It takes in a hb time series and a corresponding ScansQuality

% Set bad channels to zero (both long and short) so the ShortDistanceFilter
% will work later on. 
hb_pruned = hb;

for i=1:numel(hb_pruned)
    sChannel = find(hb_pruned(i).probe.link.ShortSeperation == 1);
   idxBadCh = find(ScansQuality(i).qMats.MeasListAct==0);
   fprintf('Scan:%i #BadChannels:%i\n',i,length(idxBadCh)/2);
   hb_pruned(i).data(:,idxBadCh) = nan; 
   hb_pruned(i).data(:,sChannel(isnan(hb_pruned(i).data(1,sChannel)))) = nan; 
end 

end
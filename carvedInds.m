function [carvedinds,H3] = carvedInds(swcfolder,RR,opt)
%CARVEDINDS Summary of this function goes here
% 
% [OUTPUTARGS] = CARVEDINDS(INPUTARGS) Explain usage here
% 
% Inputs: 
% 
% Outputs: 
% 
% Examples: 
% 
% Provide sample usage code here
% 
% See also: List related files here

% $Author: base $	$Date: 2017/02/21 16:02:02 $	$Revision: 0.1 $
% Copyright: HHMI 2017

% load locations
myswcs = dir([swcfolder,'*.swc']);
scale = 1000
for ii=1:length(myswcs)
    swcfile = fullfile(swcfolder,myswcs(ii).name);
    [swcData,offset,color, header] = loadSWC(swcfile);
    swcData(:,3:5) = swcData(:,3:5) + ones(size(swcData,1),1)*offset;
    swcData(:,3:5) = swcData(:,3:5)*scale;
    swcD{ii} = swcData;
end
%%
% find all tiles that encapsulates query data
swclocs = cat(1,swcD{:});swclocs=swclocs(:,3:5)/1e3;
pixloc = um2pix(opt,swclocs);

ed1 = unique(RR(:,1));ed1(end+1)=ed1(end)+opt.imgsiz(1);
ed2 = unique(RR(:,2));ed2(end+1)=ed2(end)+opt.imgsiz(2);
ed3 = unique(RR(:,3));ed3(end+1)=ed3(end)+opt.imgsiz(3);
H3=histcn(pixloc,ed1,ed2,ed3);
[ix,iy,iz] = ind2sub(size(H3),find(H3>0));
% H = vol3d('CData',H3>0,'texture','3D')
% alphamap('rampup');
% alphamap('rampdown'), alphamap('decrease'), alphamap('decrease')
% alphamap(.06 .* alphamap);
% find these indicies in RR

carvedinds = nan(1,length(ix));
for ii=1:length(ix)
    st1 = RR(:,1)==ed1(ix(ii));
    st2 = RR(:,2)==ed2(iy(ii));
    st3 = RR(:,3)==ed3(iz(ii));
    rrind = find(st1&st2&st3);
    carvedinds(ii) = rrind;
%     RR(rrind,:)
end

end

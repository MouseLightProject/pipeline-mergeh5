sample = '2017-09-25'
transformfile = fullfile('/nrs/mouselight/cluster/classifierOutputs/',sprintf('%s',sample),'render','transform.txt')
sample_ = sample;sample_(strfind(sample,'-'))=[];
h5config = fullfile('/groups/mousebrainmicro/home/base/CODE/MATLAB/pipeline/mergeh5/config_files/',sprintf('%s_prob%d_config_mergeh5.cfg',sample_,0));
paramsTransform = configparser(transformfile);
paramsTransform.level = paramsTransform.nl-1;
paramsH5 = configparser(h5config);
paramsH5.level = paramsH5.nl-1;
h5file = dir([paramsH5.outname,'*.h5']);
h5file = fullfile(h5file.folder,h5file.name);
%%
umloc = [70825.3, 14222.9, 35526.1];%[74907.4, 14013.1, 36651.3];%[70825.3, 14222.9, 35526.1] 
umloc = [74295.5, 13455.8, 37631.0]
%%
% find the tile that has umloc
clear ix
for ii=1:3
    ix(:,ii) = scopeloc.loc(:,ii)*1e3<umloc(ii) & (scopeloc.loc(:,ii)*1e3+imsize_um(ii))>umloc(ii);
end
[round(scopeloc.loc(find(all(ix,2)),:)*1e3) round(scopeloc.loc(find(all(ix,2)),:)*1e3+imsize_um)]
find(all(ix,2))
scopeloc.gridix(find(all(ix,2),1),:)
for ii=1:1e6;if strcmp(scopeloc.relativepaths{ii},'/2017-10-01/01/01956');break;end;end;ii
for ii=1:1e6;if strcmp(scopeloc.relativepaths{ii},'/2017-10-01/01/01983');break;end;end;ii

%%
pixloc = um2pix(paramsTransform,umloc);
% load h5 file
ROI = [1 1 1/20]*1000;
st = pixloc-ROI/2;
Io_ = permute(h5read(h5file,['/',paramsH5.h5channel],st,ROI),[2 1 3]);
Io = uint8(255*(Io_>50));
figure, imshow(max(Io_,[],3),[.1 .6]*255)

% figure,
% imshow3D(permute(Io,[2 1 3]))
%%
close all
figure(100)
cla
vol3d('CData',Io);
set(gca,'Color','k')
% daspect([1 1 3]) 
% colormap('gray')
% alphamap('rampdown'), 
% alphamap('vdown'), 
% alphamap('decrease')
% alphamap('decrease')
axis ij equal

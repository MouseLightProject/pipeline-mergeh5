function [outputArgs] = mergeh5(configfile)
%MERGEH5 Creates a single h5 file from multiple small tif/h5s. For tifs,
%provide parent folder in octree, for h5s provide envelop folder. Output is
%8 bit so make sure that threshold and parameters are set properly.
%
% [OUTPUTARGS] = MERGEH5(INPUTARGS)
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

% $Author: base $	$Date: 2016/05/16 17:25:20 $	$Revision: 0.1 $
% Copyright: HHMI 2016
if nargin<1
    configfile = './config_files/20150619_config_mergeh5.cfg'
end
addpath(genpath('./common'))
opt = configparser(configfile);

% if isfield(opt,'level')
if isfield(opt,'nl')
    opt.level=opt.nl-1;
else
    error('No level was found')
end
%%
% set block size to a fraction of image size to maximize speed
[~,~,fileext] = fileparts(opt.ext);
if strcmp(fileext,'.h5')
    myh5 = dir(fullfile(opt.inputfolder,'*.h5'));
    info = h5info(fullfile(opt.inputfolder,myh5(1).name));
    imgsiz = info.Datasets.Dataspace.Size;
    opt.imgsiz = imgsiz;
    blocksize = imgsiz/2;
    while any(blocksize>[128 128 128])
        div = blocksize>[128 128 128];
        blocksize = blocksize./(div+1);
    end
    outsiz = opt.imgsiz*2^(opt.level);
else
    mytif = dir(fullfile(opt.inputfolder,'*.tif'));
    info = imfinfo(fullfile(opt.inputfolder,mytif(1).name), 'tif');
    imgsiz = double([info(1).Width info(1).Height length(info)]);
    opt.imgsiz = imgsiz;
    blocksize = imgsiz/2;
    while any(blocksize>[128 128 128])
        div = blocksize>[128 128 128];
        blocksize = blocksize./(div+1);
    end
    outsiz = opt.imgsiz*2^(opt.level);
end
% if any(opt.imgsiz~=whd)
%     opt.imgsiz = whd;
%     warning(sprintf('tile size set to [%d %d %d]',opt.imgsiz))
% end
%%
% get sequence
args.level = opt.level;
args.ext = opt.ext;
if exist(opt.seqtemp, 'file') == 2
    % load file directly
else
    args.fid = fopen(opt.seqtemp,'w');
    recdir(opt.inputfolder,args)
end
fid=fopen(opt.seqtemp,'r');
myfiles = textscan(fid,'%s');
myfiles = myfiles{1};
fclose(fid)
%%
RR = zeros(size(myfiles,1),6);
SS = zeros(size(myfiles,1),3);
parfor idx=1:size(myfiles,1)
    %%
    tmps = strsplit(myfiles{idx},filesep);
    seq = [tmps{end-opt.level:end-1}];
    idxtile = str2num(seq);
    lenseq = ceil(log10(idxtile));
    st = [0 0 0];
    bin = opt.imgsiz'*2.^[lenseq-1:-1:0];
    indxyz = zeros(lenseq,3);
    for iseq = 1:lenseq
        is = seq(iseq);
        temp = fliplr(dec2bin(str2num(is)-1,3)); % coordinates wrto xyz
        for ii=1:3 % x y z
            indxyz(iseq,ii) = str2num(temp(ii));
            st(ii) = st(ii) + bin(ii,iseq)*indxyz(iseq,ii);
        end
    end
    RR(idx,:)=[st st+opt.imgsiz];
    SS(idx,:) = [rem(floor(idxtile/10^(opt.level-1)),10) rem(floor(idxtile/10^(opt.level-2)),10) rem(floor(idxtile/10^(opt.level-3)),10)];
end
%%
myh5prob =sprintf('/%s',opt.h5channel);
myouth5 = sprintf('%s_lev-%d.h5',opt.outname,opt.level);
%%
clc
% myouth5 = 'mytest.h5';
% unix(sprintf('rm %s',myouth5))
h5create(myouth5,myh5prob,outsiz,'Datatype','single','ChunkSize',blocksize,'Deflate',1)
h5create(myouth5,[myh5prob,'_props/origin'], [1 3]);
h5create(myouth5,[myh5prob,'_props/spacing'], [1 3]);
h5create(myouth5,[myh5prob,'_props/level'], [1]);
h5create(myouth5,[myh5prob,'_props/ROI'], size(RR));
h5create(myouth5,[myh5prob,'_props/telapsed'], size(RR,1));
%%
h5write(myouth5,[myh5prob,'_props/origin'], [opt.ox opt.oy opt.oz]);
h5write(myouth5,[myh5prob,'_props/spacing'], [opt.sx opt.sy opt.sz]);
h5write(myouth5,[myh5prob,'_props/level'], [opt.nl-1]);
h5write(myouth5,[myh5prob,'_props/ROI'], RR);
%%
telapsed = zeros(size(myfiles,1),1);
tstart = tic;
for idx=1:size(myfiles)
    %% write into big file
    st = RR(idx,1:3);
    if strcmp(fileext,'.h5')
        data = single(h5read(myfiles{idx},['/',info.Datasets.Name]));
        data = uint8(((data+1).*(single(data>opt.maskThr)/256))-1);
    else
        data = permute(single(deployedtiffread(myfiles{idx})),[2 1 3]);
        data = uint8(((data+1).*(single(data>opt.maskThr)/256))-1);
    end
    h5write(myouth5,myh5prob,data,st+[1 1 1],size(data),[1 1 1])
    ttoc = round(toc(tstart));
    telapsed(idx) = ttoc;
    disp(sprintf('Elapsed time: %dsecs / Working on %d out of %d / St: [%d %d %d]',ttoc,idx,size(myfiles,1),st(1),st(2),st(3)))
end
h5write(myouth5,[myh5prob,'_props/telapsed'], telapsed);

%%
if opt.viz
    figure,
    hold on
    plotyy(1:length(telapsed),telapsed,1:length(telapsed)-1,diff(telapsed))
    title(num2str(max(telapsed)))
    legend('total time','iter time')
    % load seq
    BB = RR(:,[1 4 2 5 3 6])+1;
    %%
end

end

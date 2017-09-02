function mergeh5_chunks(configfile,iiC,jjC,kkC,numchunk)
%MERGEH5 Creates multiple h5 files from multiple small tif/h5s. For tifs,
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
setmask=0;
if nargin<1
    configfile = './config_files/20170810_prob1_config_mergeh5.cfg'
elseif nargin==1
    
elseif nargin ==5
    iiC = str2double(iiC);
    jjC = str2double(jjC);
    kkC = str2double(kkC);
    numchunk = str2double(numchunk);
end
if ~isdeployed
    addpath(genpath('./common'))
end
opt = configparser(configfile);

% if isfield(opt,'level')
if isfield(opt,'nl')
    opt.level=opt.nl-1;
else
    error('No level was found')
end
if isfield(opt,'numCPU')
    numCPU = opt.numCPU;
else
    numCPU = feature('numcores');
end
poolobj = gcp('nocreate'); % If no pool, do not create new one.
if isempty(poolobj)
    poolsize = 0;
    parpool(numCPU)
else
    poolsize = poolobj.NumWorkers
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
    mkdir(fileparts(opt.seqtemp))
    args.fid = fopen(opt.seqtemp,'w');
    recdir(opt.inputfolder,args)
end
fid=fopen(opt.seqtemp,'r');
myfiles = textscan(fid,'%s');
myfiles = myfiles{1};
fclose(fid);
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
%% mask a region based on an given swc
if setmask
    swcfolder = '/groups/mousebrainmicro/home/base/CODE/MATLAB/pipeline/mergeh5/sampleswcs/'
    [carvedinds,H3] = carvedInds(swcfolder,RR,opt);
    myfiles = myfiles(carvedinds);
    RR = RR(carvedinds,:);
    SS = SS(carvedinds,:);
end
%%
myh5prob =sprintf('/%s',opt.h5channel);
%% go 3 level deep
if 0
    numchunk = 8;
    for ii=1:numchunk
        for jj=1:numchunk
            for kk=1:numchunk
                %%
                myouth5 = sprintf('%s_lev-%d_chunk-%d%d%d_%d%d%d.h5',opt.outname,opt.level,ii,jj,kk,numchunk,numchunk,numchunk);
                %myouth5 = sprintf('/data/lev-%d_chunk-%d%d%d_%d%d%d.h5',opt.level,ii,jj,kk,numchunk,numchunk,numchunk);
                mkdir(fileparts(myouth5))
                h5create(myouth5,myh5prob,outsiz,'Datatype','single','ChunkSize',blocksize,'Deflate',1)
                h5create(myouth5,[myh5prob,'_props/origin'], [1 3]);
                h5create(myouth5,[myh5prob,'_props/spacing'], [1 3]);
                h5create(myouth5,[myh5prob,'_props/level'], [1]);
                %h5create(myouth5,[myh5prob,'_props/ROI'], size(RR));
                
                h5create(myouth5,[myh5prob,'_props/telapsed'], size(RR,1));
                h5write(myouth5,[myh5prob,'_props/origin'], [opt.ox opt.oy opt.oz]);
                h5write(myouth5,[myh5prob,'_props/spacing'], [opt.sx opt.sy opt.sz]);
                h5write(myouth5,[myh5prob,'_props/level'], [opt.nl-1]);
                % h5write(myouth5,[myh5prob,'_props/ROI'], RR);, only write
                % where the data is copied
                %%
                idxiijjkk = find(SS(:,1)==ii&SS(:,2)==jj&SS(:,3)==kk);
                if isempty(idxiijjkk)
                    continue
                end
                h5create(myouth5,[myh5prob,'_props/ROI'], size(RR(idxiijjkk,:)));
                h5write(myouth5,[myh5prob,'_props/ROI'], RR(idxiijjkk,:));
                h5create(myouth5,[myh5prob,'_props/ROISS'], size(SS(idxiijjkk,:)));
                h5write(myouth5,[myh5prob,'_props/ROISS'], SS(idxiijjkk,:));
                % read all
                %% parread()
                theseinds = idxiijjkk(:)';
                Itempsub = cell(1,length(theseinds));
                parfor idx = 1:length(theseinds)
                    if strcmp(fileext,'.h5')
                        data = single(h5read(myfiles{theseinds(idx)},['/',info.Datasets.Name]));
                        data = uint8(((data+1).*(single(data>opt.maskThr)/256))-1);
                    else
                        data = permute(single(deployedtiffread(myfiles{theseinds(idx)})),[2 1 3]);
                        data = uint8(((data+1).*(single(data>opt.maskThr)/256))-1);
                    end
                    %st = RR(theseinds(idx),1:3)-bbox(1:3);
                    Itempsub{idx} = data;
                end
                %%
                telapsed = zeros(1,length(theseinds));
                for idx = 1:length(theseinds)
                    % write into big file
                    st = RR(theseinds(idx),1:3);
                    data = Itempsub{idx};
                    tstart = tic;
                    h5write(myouth5,myh5prob,data,st+[1 1 1],size(data),[1 1 1])
                    ttoc = round(toc(tstart));
                    telapsed(idx) = ttoc;
                    %         disp(sprintf('Elapsed time: %dsecs / Working on %d out of %d / St: [%d %d %d]',ttoc,idx,size(myfiles,1),st(1),st(2),st(3)))
                end
                %%
                if 0
                    %%
                    rr = RR(idxiijjkk,:);
                    bbox = [min(rr(:,1:3)) max(rr(:,4:6))]
                    tmpsize = bbox(4:6)-bbox(1:3);
                    Itemp = zeros(tmpsize(1),tmpsize(2),tmpsize(3),'uint8');
                    for idx = 1:length(theseinds)
                        st = RR(theseinds(idx),1:3)-bbox(1:3);
                        Itemp(st(1)+1:st(1)+opt.imgsiz(1),...
                            st(2)+1:st(2)+opt.imgsiz(2),...
                            st(3)+1:st(3)+opt.imgsiz(3)) = Itempsub{idx};
                    end
                    %%
                    myh5prob =sprintf('/%s',opt.h5channel);
                    myouth5 = sprintf('/data/chunk.h5');
                    
                    h5create(myouth5,myh5prob,outsiz,'Datatype','single','ChunkSize',blocksize,'Deflate',1)
                    tic
                    h5write(myouth5,myh5prob,Itemp,bbox(1:3)+[1 1 1],size(Itemp),[1 1 1])
                    sprintf('write chunk in %f',toc)
                    
                end
                %%
                
            end
        end
    end
else
    myouth5 = sprintf('%s_lev-%d_chunk-%d%d%d_%d%d%d_masked-%d.h5',opt.outname,opt.level,1,1,1,1,1,1,setmask);
    
    %myouth5 = sprintf('/data/lev-%d_chunk-%d%d%d_%d%d%d.h5',opt.level,ii,jj,kk,numchunk,numchunk,numchunk);
    mkdir(fileparts(myouth5))
    h5create(myouth5,myh5prob,outsiz,'Datatype','single','ChunkSize',blocksize,'Deflate',1)
    h5create(myouth5,[myh5prob,'_props/origin'], [1 3]);
    h5create(myouth5,[myh5prob,'_props/spacing'], [1 3]);
    h5create(myouth5,[myh5prob,'_props/level'], [1]);
    %h5create(myouth5,[myh5prob,'_props/ROI'], size(RR));
    
    h5create(myouth5,[myh5prob,'_props/telapsed'], size(RR,1));
    h5write(myouth5,[myh5prob,'_props/origin'], [opt.ox opt.oy opt.oz]);
    h5write(myouth5,[myh5prob,'_props/spacing'], [opt.sx opt.sy opt.sz]);
    h5write(myouth5,[myh5prob,'_props/level'], [opt.nl-1]);
    % h5write(myouth5,[myh5prob,'_props/ROI'], RR);, only write
    % where the data is copied
    %%
    % idxiijjkk = find(SS(:,1)==iiC&SS(:,2)==jjC&SS(:,3)==kkC);
    idxiijjkk=1:length(myfiles);
    if isempty(idxiijjkk)
        return
    end
    h5create(myouth5,[myh5prob,'_props/ROI'], size(RR(idxiijjkk,:)));
    h5write(myouth5,[myh5prob,'_props/ROI'], RR(idxiijjkk,:));
    % read all
    
    %% parread()
    theseinds = idxiijjkk(:)';
    telapsed = zeros(1,length(theseinds));
    %Itempsub = cell(1,length(theseinds));
    Itempsub = cell(1,numCPU);
    
    iters = 0:numCPU:length(theseinds);
    kk=length(iters)-1;
    for ii=1:kk
        sprintf('%d out of %d',ii,kk)
        ticparread = tic;
        % read in paralel
        parfor idx = 1:numCPU%length(theseinds)
            if strcmp(fileext,'.h5')
                data = single(h5read(myfiles{iters(ii)+idx},['/',info.Datasets.Name]));
                data = uint8(((data+1).*(single(data>opt.maskThr)/256))-1);
            else
                data = permute(single(deployedtiffread(myfiles{theseinds(idx)})),[2 1 3]);
                data = uint8(((data+1).*(single(data>opt.maskThr)/256))-1);
            end
            %st = RR(theseinds(idx),1:3)-bbox(1:3);
            Itempsub{idx} = data;
        end
        elapseread = toc(ticparread);
        %%
        %%% write sequential
        for idx = 1:numCPU
            % write into big file
            st = RR(iters(ii)+idx,1:3);
            data = Itempsub{idx};
            tstart = tic;
            h5write(myouth5,myh5prob,data,st+[1 1 1],size(data),[1 1 1])
            ttoc = round(toc(tstart));
            telapsed(iters(ii)+idx) = ttoc;
        end
        sprintf('block idx: %d / %d, R: %d, W: %d, maxW: %f',ii,kk,round(elapseread),...
            round(sum(telapsed(iters(ii)+1:iters(ii)+numCPU))),max(telapsed(iters(ii)+1:iters(ii)+numCPU)))

    end
    %%
    for idx = iters(end)+1: length(theseinds)
        % write into big file
        st = RR(idx,1:3);
        if strcmp(fileext,'.h5')
            data = single(h5read(myfiles{idx},['/',info.Datasets.Name]));
            data = uint8(((data+1).*(single(data>opt.maskThr)/256))-1);
        else
            data = permute(single(deployedtiffread(myfiles{theseinds(idx)})),[2 1 3]);
            data = uint8(((data+1).*(single(data>opt.maskThr)/256))-1);
        end
        tstart = tic;
        h5write(myouth5,myh5prob,data,st+[1 1 1],size(data),[1 1 1])
        ttoc = round(toc(tstart));
        telapsed(idx) = ttoc;
    end
end
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
function deployment
% mcc -m -R -nojvm -v mergeh5_chunks.m -d /groups/mousebrainmicro/home/base/CODE/MATLAB/compiledfunctions/mergeh5_chunks -a ./common
%
mergeh5_chunks('/groups/mousebrainmicro/home/base/CODE/MATLAB/recontree/config_files/20150619_config_mergeh5.cfg','1','1','1','8')

%%
brain = '2015-06-19';
configfile = '/groups/mousebrainmicro/home/base/CODE/MATLAB/recontree/config_files/20161025_config_mergeh5_ch1.cfg'
configfile = '/groups/mousebrainmicro/home/base/CODE/MATLAB/recontree/config_files/20150619_config_mergeh5.cfg'
opt = configparser(configfile);
if isfield(opt,'nl')
    opt.level=opt.nl-1;
else
    error('No level was found')
end
% get sequence
args.level = opt.level;
args.ext = opt.ext;
if exist(opt.seqtemp, 'file') == 2
    % load file directly
else
    mkdir(fileparts(opt.seqtemp))
    args.fid = fopen(opt.seqtemp,'w');
    recdir(opt.inputfolder,args)
end
fid=fopen(opt.seqtemp,'r');
myfiles = textscan(fid,'%s');
myfiles = myfiles{1};
fclose(fid);
%%
% mergeh5_chunks(configfile,iiC,jjC,kkC,numchunk)
% mergeh5_chunks(configfile,'2','1','1','8')
% outputfold = '/groups/mousebrainmicro/mousebrainmicro/cluster/Stitching/160718/Descriptors/'
numcores = 4;
myfile = sprintf('mergeh5chunks_%s_%s_ch0.sh',brain,date)
compiledfunc = '/groups/mousebrainmicro/home/base/CODE/MATLAB/compiledfunctions/mergeh5_chunks/mergeh5_chunks'

%find number of random characters to choose from
s = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789';
numRands = length(s);
%specify length of random string to generate
sLength = 10;
%-o /dev/null
esttime = 30*60;

numchunk = 8;
[subs1 subs2 subs3] = ndgrid(1:8,1:8,1:8);
subs = [subs3(:) subs2(:) subs1(:)];
% if 1
%     myfiles = dir(fullfile(matfolder,'pointmatches','*.mat'));
%     doneinds = cellfun(@(x) str2num(x(1:5)),{myfiles.name});
%     [finished,bb] = min(pdist2((inds+1)',doneinds(:)),[],2);finished = ~finished;
%     % [finished,bb] = min(pdist2((inds+1)',find(cellfun(@isempty,regpts))'),[],2)
% else
%     finished = zeros(1,length(inds)-1);
% %     finished(402:494) = 1;
% %     finished = ~finished
% end
%%
% mergeh5_chunks(configfile,'1','2','3','8')
%%
fid = fopen(myfile,'w');
for ii=1:size(subs,1)
    %%
    %generate random string
    randString = s( ceil(rand(1,sLength)*numRands) );
    name = sprintf('mh5_%05d-%s',ii,randString);
    outargs = sprintf('''%s %s %d %d %d %d> output.log''',compiledfunc,configfile,subs(ii,1),subs(ii,2),subs(ii,3),numchunk);
    mysub = sprintf('qsub -pe batch %d -l d_rt=%d -N %s -j y -o ~/logs -b y -cwd -V %s\n',numcores,esttime,name,outargs);
    % check if output files exists
    myouth5 = sprintf('%s_lev-%d_chunk-%d%d%d_%d%d%d.h5',opt.outname,opt.level,subs(ii,1),subs(ii,2),subs(ii,3),numchunk,numchunk,numchunk);
    if exist(myouth5,'file')
        continue
    else
        fwrite(fid,mysub);
    end
end
unix(sprintf('chmod +x %s',myfile));

end

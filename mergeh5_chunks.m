function mergeh5_chunks(configuration_file_name)
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

    % Read the configuration file, break out the things in it
    opt = configparser(configuration_file_name);
    input_folder_name = opt.inputfolder ;
    %nl = opt.nl ;
    if isfield(opt,'numCPU')
        core_count_requested = opt.numCPU;
    else
        core_count_requested = feature('numcores');
    end
    output_file_name = opt.outname ;
    input_file_names_list_file_path = opt.seqtemp ;
    h5_dataset_name = opt.h5channel ;
    input_file_name_ending = opt.ext ;
    mask_threshold = opt.maskThr ;
    do_visualize = opt.viz ;

    % inputfolder = '/nrs/mouselight/Users/mluser/2018-05-23-prob'
    % # output h5 name
    % outname = '/scratch/classifierOutputs/2018-04-13/20180413_prob0/20180413_prob0'
    % # filelist sequence
    % seqtemp = '/scratch/classifierOutputs/2018-04-13/20180413_prob0/20180413_prob0-seq0.txt'
    % # copy scratch to /nrs/mouselight/cluster/classifierOutputs/2018-04-13/20180413_prob0

    optTransform = configparser(fullfile(input_folder_name, 'transform.txt'));
    ox = optTransform.ox ;
    oy = optTransform.oy ;
    oz = optTransform.oz ;
    sx = optTransform.sx ;
    sy = optTransform.sy ;
    sz = optTransform.sz ;
    nl = optTransform.nl ;
    
%     % for thesefields = {'ox','oy','oz','sx','sy','sz','nl'}
%     %     (thesefields{1}) = optTransform.(thesefields{1});
%     % end

    level = nl-1 ;
    poolobj = gcp('nocreate'); % If no pool, do not create new one.
    if isempty(poolobj)
        %poolsize = 0;
        parpool([1 core_count_requested])
    else
        %poolsize = poolobj.NumWorkers
    end
    %%
    % set block size to a fraction of image size to maximize speed
    [~,~,fileext] = fileparts(input_file_name_ending);
    % leaf image stack size
    if strcmp(fileext,'.h5')
        myh5 = dir(fullfile(input_folder_name,'*.h5'));
        info = h5info(fullfile(input_folder_name,myh5(1).name));
        imgsiz = info.Datasets.Dataspace.Size;
        if length(imgsiz)>3 ,
            imgsiz = imgsiz(end-2:end) ;  
        end
    else
        mytif = dir(fullfile(input_folder_name,'default.0.tif'));
        info = imfinfo(fullfile(input_folder_name,mytif(1).name), 'tif');
        imgsiz = double([info(1).Width info(1).Height length(info)]);
    end

    %imgsiz = imgsiz;
    blocksize = imgsiz/2;
    while any(blocksize>[128 128 128])
        div = blocksize>[128 128 128];
        blocksize = blocksize./(div+1);
    end
    outsiz = imgsiz*2^(level);

    %%
    % get sequence
    args = struct() ;
    args.level = level;
    args.ext = input_file_name_ending;
    if ~exist(input_file_names_list_file_path, 'file') ,
        fprintf('Creating file: %s\n',input_file_names_list_file_path) ;
        parent_folder_path = fileparts(input_file_names_list_file_path) ;
        if ~exist(parent_folder_path, 'file') ,
            mkdir(parent_folder_path) ;
        end
        args.fid = fopen(input_file_names_list_file_path, 'w') ;
        recdir(input_folder_name,args)
    end
    fid=fopen(input_file_names_list_file_path,'r');
    myfiles_raw = textscan(fid,'%s');
    myfiles = myfiles_raw{1};
    fclose(fid);
    %%
    RR = zeros(size(myfiles,1),6);
    %SS = zeros(size(myfiles,1),3);
    parfor idx=1:size(myfiles,1)
        %%
        tmps = strsplit(myfiles{idx},filesep);
        seq = [tmps{end-level:end-1}];
        idxtile = str2num(seq);  %#ok<ST2NM>
        lenseq = ceil(log10(idxtile));
        st = [0 0 0];
        bin = imgsiz'*2.^(lenseq-1:-1:0);
        indxyz = zeros(lenseq,3);
        for iseq = 1:lenseq
            is = seq(iseq);
            temp = fliplr(dec2bin(str2num(is)-1,3)); %#ok<ST2NM>  % coordinates wrto xyz
            for ii=1:3 % x y z
                indxyz(iseq,ii) = str2num(temp(ii));  %#ok<ST2NM>
                st(ii) = st(ii) + bin(ii,iseq)*indxyz(iseq,ii);
            end
        end
        RR(idx,:)=[st st+imgsiz];
        %SS(idx,:) = [rem(floor(idxtile/10^(level-1)),10) rem(floor(idxtile/10^(level-2)),10) rem(floor(idxtile/10^(level-3)),10)];
    end
    %% mask a region based on an given swc
    % if 0
    %     swcfolder = '/groups/mousebrainmicro/home/base/CODE/MATLAB/pipeline/mergeh5/sampleswcs/'
    %     [carvedinds,H3] = carvedInds(swcfolder,RR,opt);
    %     myfiles = myfiles(carvedinds);
    %     RR = RR(carvedinds,:);
    %     SS = SS(carvedinds,:);
    % end
    %%
    myh5prob =sprintf('/%s',h5_dataset_name);

    %%
    %output_file_name = sprintf('%s_lev-%d_chunk-%d%d%d_%d%d%d_masked-%d.h5',output_folder_name,level,1,1,1,1,1,1,0);

    %myouth5 = sprintf('/data/lev-%d_chunk-%d%d%d_%d%d%d.h5',level,ii,jj,kk,numchunk,numchunk,numchunk);
    myouth5_parent_folder_path = fileparts(output_file_name) ;
    if ~exist(myouth5_parent_folder_path, 'file') ,
        mkdir(fileparts(output_file_name))
    end
    if exist(output_file_name, 'file') ,
        error('Target file %s already exists', output_file_name) ;
    end
    h5create(output_file_name,myh5prob,outsiz,'Datatype','single','ChunkSize',blocksize,'Deflate',3)
    h5create(output_file_name,[myh5prob,'_props/origin'], [1 3]);
    h5create(output_file_name,[myh5prob,'_props/spacing'], [1 3]);
    h5create(output_file_name,[myh5prob,'_props/level'], 1);
    %h5create(myouth5,[myh5prob,'_props/ROI'], size(RR));

    h5create(output_file_name,[myh5prob,'_props/telapsed'], size(RR,1));
    h5write(output_file_name,[myh5prob,'_props/origin'], [ox oy oz]);
    h5write(output_file_name,[myh5prob,'_props/spacing'], [sx sy sz]);
    h5write(output_file_name,[myh5prob,'_props/level'], nl-1);
    % h5write(myouth5,[myh5prob,'_props/ROI'], RR);, only write
    % where the data is copied
    %%
    % idxiijjkk = find(SS(:,1)==iiC&SS(:,2)==jjC&SS(:,3)==kkC);
    idxiijjkk=1:length(myfiles);
    if isempty(idxiijjkk)
        return
    end
    h5create(output_file_name,[myh5prob,'_props/ROI'], size(RR(idxiijjkk,:)));
    h5write(output_file_name,[myh5prob,'_props/ROI'], RR(idxiijjkk,:));
    % read all

    %% parread()
    theseinds = idxiijjkk(:)';
    telapsed = zeros(1,length(theseinds));
    %Itempsub = cell(1,length(theseinds));
    numBatch = 4*core_count_requested;
    Itempsub = cell(1,numBatch);

    iters = 0:numBatch:length(theseinds) ;
    kk=length(iters)-1;
    for ii=1:kk
        %%
        sprintf('%d out of %d',ii,kk)
        ticparread = tic;
        % read in paralel
        parfor idx = 1:numBatch%length(theseinds)
            if strcmp(fileext,'.h5')
                data = single(h5read(myfiles{iters(ii)+idx},['/',info.Datasets.Name]));  %#ok<PFBNS>
                data = uint8(((data+1).*(single(data>mask_threshold)/256))-1);
                if ndims(data)>3 ,
                    original_size = size(data) ;
                    new_size = original_size(end-2:end) ;
                    data = reshape(data, new_size) ;  
                end
            else
                data = permute(single(deployedtiffread(myfiles{iters(ii)+idx})),[2 1 3]);
                data = (data+1).*single(data>mask_threshold);
                data(~data) = mask_threshold;
                data = data/256-1;
                data = uint8(data);
            end
            %st = RR(theseinds(idx),1:3)-bbox(1:3);
            Itempsub{idx} = data;
        end
        elapseread = toc(ticparread);
        %%
        %%% write sequential
        for idx = 1:numBatch
            % write into big file
            st = RR(iters(ii)+idx,1:3);
            data = Itempsub{idx};
            tstart = tic;
            h5write(output_file_name,myh5prob,data,st+[1 1 1],size(data),[1 1 1])
            ttoc = (toc(tstart));
            telapsed(iters(ii)+idx) = ttoc;
        end

        sprintf('block idx: %d / %d, R: %d, W: %d, maxW: %f',ii,kk,round(elapseread),...
            round(sum(telapsed(iters(ii)+1:iters(ii)+numBatch))),max(telapsed(iters(ii)+1:iters(ii)+numBatch)))
    end

    for idx = iters(end)+1: length(theseinds)
        % write into big file
        st = RR(idx,1:3);
        if strcmp(fileext,'.h5')
            data = single(h5read(myfiles{idx},['/',info.Datasets.Name]));
            data = uint8(((data+1).*(single(data>mask_threshold)/256))-1);
            if ndims(data)>3 ,
                original_size = size(data) ;
                new_size = original_size(end-2:end) ;
                data = reshape(data, new_size) ;
            end                
        else
            data = permute(single(deployedtiffread(myfiles{theseinds(idx)})),[2 1 3]);
            data = (data+1).*single(data>mask_threshold);
            data(~data) = mask_threshold;
            data = data/256-1;
            data = uint8(data);
        end
        tstart = tic;
        h5write(output_file_name,myh5prob,data,st+[1 1 1],size(data),[1 1 1])
        ttoc = round(toc(tstart));
        telapsed(idx) = ttoc;
    end
    % copy output to target location

    % delete the scratch location

    %%
    if do_visualize
        figure,
        hold on
        plotyy(1:length(telapsed),telapsed,1:length(telapsed)-1,diff(telapsed)) %#ok<PLOTYY>
        title(num2str(max(telapsed)))
        legend('total time','iter time')
        % load seq
        %BB = RR(:,[1 4 2 5 3 6])+1;
        %%
    end

end  % main function



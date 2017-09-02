configfile = './config_files/20150619_config_mergeh5.cfg'
opt = configparser(configfile);
opt.level=opt.nl-1;
args.level = opt.level;
%%
args.ext = 'tif';
args.level=4;
opt.seqtemp = '/nrs/mouselight/2015-06-19-johan-full/filelist-2.txt'
opt.inputfolder = '/nrs/mouselight/2015-06-19-johan-full/2'
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

%%
s = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789';
compiledfunc = '/groups/mousebrainmicro/home/base/CODE/syglass/syCacher'
%find number of random characters to choose from
numRands = length(s);
%specify length of random string to generate
sLength = 10;
numcores = 4;
timelim = 5*60
fold=1
mysh = sprintf('syglassrun-%d.sh',fold);

fid = fopen(mysh,'w');
for ix=1:8
    % syGlass -i <path to the source folder> -o <path to the destination folder> -t -x 10000
    for idx = 1:64
        %     sprintf('syGlass -i <path to the source folder> -o <path to the destination folder> -t -x 10000')
        
        [bb,aa]=ind2sub([8 8],idx);
        randString = s( ceil(rand(1,sLength)*numRands) );
        
        infold = sprintf('/nrs/mouselight/2015-06-19-johan-full/%d/%d/%d/%d',fold,ix,aa,bb);
        outfold = sprintf('/nrs/mouselight/2015-06-19-johan-full-syg/%d/%d/%d/%d/',fold,ix,aa,bb);
        %     mkdir(outfold)
        if ~exist(fullfile(infold,sprintf('default.0.tif')))
            
        else
            if exist(fullfile(outfold,sprintf('default.0.tif')))
                [ix idx]
            else
                
                %     outfile `= fullfile(outfolder,sprintf('%s_idx-%05d_stxyzendxyz-%d_%d_%d_%d_%d_%d.txt',bb,idx,BB(1:2:end),BB(2:2:end)));
                name = sprintf('cache_%05d-%s',(ix-1)*8+idx,randString);
                args = sprintf('''%s -i %s -o %s -t -x 10000''',compiledfunc,infold,outfold);
                mysub = sprintf('LD_LIBRARY_PATH=/misc/local/gcc-5.1.0/lib64/ qsub -pe batch %d -l d_rt=%d -N %s -j y -o /dev/null -b y -cwd -V %s\n',numcores,timelim,name,args);
                fwrite(fid,mysub);
            end
        end
    end
end
%%
% bsub -n1 -We 1 -J t-00001 -o /dev/null 'python singleThreadedCacher.py  /groups/turaga/home/moreheadm /groups/turaga/home/moreheadm/outtest/'
% find -type d -printf '%d\t%P\n' | sort -r -nk4 | cut -f2- > listfolders
foldind = 2
% infold = sprintf('/nrs/mouselight/SAMPLES/2017-06-28-2/%d',foldind);
infold = '/nrs/mouselight/SAMPLES/2017-06-28-2/';
args = '-type d -printf ''%d\t%P\n'' | sort -r -nk1 | cut -f2- >';
outfile = fullfile(infold,'listfolders');
unix(sprintf('find %s %s %s',infold,args,outfile))
%%
fid=fopen(outfile,'r');
myfiles = textscan(fid,'%s');
myfiles = myfiles{1};
fclose(fid)
%%
samp = '2017-06-28-2-sy'
mysh = sprintf('syglassrun-%s-%d.sh',samp,foldind);

outfold = fullfile('/nrs/mouselight/SAMPLES/',samp);
fid = fopen(mysh,'w');
%%
for ii=351235:length(myfiles)
    infile = fullfile(infold,myfiles{ii});
%     if all(myfiles{ii}(1:3)=='ktx')
%         continue
%     end
    outfile = fullfile(outfold,myfiles{ii});
    logfile = sprintf('%s/%d.txt',fullfile('/groups/mousebrainmicro/mousebrainmicro/LOG/',samp),ii);
    
    mkdir(outfile)
    unix(sprintf('chmod g+rwx %s',outfile));
%     myarg = sprintf('bsub -n1 -We 1 -J t-%05d -o /dev/null ''python /groups/mousebrainmicro/mousebrainmicro/Software/syGlassConverter/singleThreadedCacher.py  %s %s %s''\n',ii,infile,outfile,logfile);
    myarg = sprintf('bsub -n1 -We 1 -J t-%05d -o /dev/null ''python /groups/mousebrainmicro/mousebrainmicro/Software/syGlassConverter/singleThreadedCacher.py  %s %s''\n',ii,infile,outfile);
    fwrite(fid,myarg);
    
end
fclose(fid)
%%



s = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789';
compiledfunc = '/groups/mousebrainmicro/home/base/CODE/syglass/syCacher'
%find number of random characters to choose from
numRands = length(s);
%specify length of random string to generate
sLength = 10;
numcores = 4;
timelim = 5*60
fold=1
mysh = sprintf('syglassrun-%d.sh',fold);

fid = fopen(mysh,'w');
for ix=1:8
    % syGlass -i <path to the source folder> -o <path to the destination folder> -t -x 10000
    for idx = 1:64
        %     sprintf('syGlass -i <path to the source folder> -o <path to the destination folder> -t -x 10000')
        
        [bb,aa]=ind2sub([8 8],idx);
        randString = s( ceil(rand(1,sLength)*numRands) );
        
        infold = sprintf('/nrs/mouselight/SAMPLES/2017-06-28-2/%d/%d/%d/%d',fold,ix,aa,bb);
        outfold = sprintf('/nrs/mouselight/2015-06-19-johan-full-syg/%d/%d/%d/%d/',fold,ix,aa,bb);
        %     mkdir(outfold)
        if ~exist(fullfile(infold,sprintf('default.0.tif')))
            
        else
            if exist(fullfile(outfold,sprintf('default.0.tif')))
                [ix idx]
            else
                
%                 
%                 %     outfile = fullfile(outfolder,sprintf('%s_idx-%05d_stxyzendxyz-%d_%d_%d_%d_%d_%d.txt',bb,idx,BB(1:2:end),BB(2:2:end)));
%                 name = sprintf('cache_%05d-%s',(ix-1)*8+idx,randString);
%                 args = sprintf('''%s -i %s -o %s -t -x 10000''',compiledfunc,infold,outfold);
%                 mysub = sprintf('LD_LIBRARY_PATH=/misc/local/gcc-5.1.0/lib64/ qsub -pe batch %d -l d_rt=%d -N %s -j y -o /dev/null -b y -cwd -V %s\n',numcores,timelim,name,args);
%                 fwrite(fid,mysub);
            end
        end
    end
end



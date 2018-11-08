# Usage: 
mergeh5_chunks(configfile)  

# Inputs:  
**configfile:**  txt file that has all the parameters used in conversion  
## Fields:  
**inputfolder:** input render folder that has the h5 octree  
**outname:** fullpath of output h5  
**seqtemp:** temporary file used to cache input tile paths that are found recursively in inputfolder at a certain depth. Depth is set to leaf node by default. useful to rerun conversion with different setting.   
**h5channel:** color channel used in wrapper  
**ext:** post-tag for naming each channel. we used 0/1 convention for our datasets.  
**maskThr:** threshold to make data sparse, useful with compression  
**viz:** 0 by default. for debugging purposes.   
**transform.txt arguments:** [ox,oy,ox,sz,sy,sz,nl]: copied from transformed.txt. @TODO@ will get rid of this in the next version as we already pass inputfolder.  
**numCPU:** number of cores used for paralel read. Useful to set a minumum number. Internaly we check max hardware cores, i.e. min(numCPU,maximum available CPU)  
  
Sample config file:
```
# parameters

# render output to be merged
inputfolder = '/nrs/mouselight/SAMPLES/2018-10-01-prob'

# output h5 name
outname = '/scratch/classifierOutputs/2018-10-01/20181001_prob0/20181001_prob0'
# filelist sequence
seqtemp = '/scratch/classifierOutputs/2018-10-01/20181001_prob0/20181001_prob0-seq0.txt'
# copy scratch to /nrs/mouselight/cluster/classifierOutputs/2018-08-01/20180801_prob0

h5channel = 'prob0'
ext = '0.h5'

maskThr = 6500
viz = 0;

# from transform.txt of the render file
ox: 69445196
oy: 12917408
oz: 30199640
sx: 16173.689102564102
sy: 16073.159226190477
sz: 64319.413043478264
nl: 7

numCPU = 40
```

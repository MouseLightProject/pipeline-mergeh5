function RR = RR_from_laden_chunk_file_paths(laden_chunk_file_paths, zoom_step_count, chunk_size_ijk)
    laden_chunk_count = length(laden_chunk_file_paths) ;  % "laden" meaning non-empty, I.e. the chunks actually represented in the 'octree'
%     RR = zeros(laden_chunk_count, 6) ;  % one row for each (actually extant) chunk at the highest zoom level
%     for laden_chunk_index = 1:laden_chunk_count ,
%         tmps = strsplit(laden_chunk_file_paths{laden_chunk_index},filesep);
%         seq = [tmps{end-zoom_step_count:end-1}];  % octree path as an array of chars, each char a digit on 1...8
%         idxtile = str2num(seq);  %#ok<ST2NM>  % octree path as an interger, the path is the digit sequence
%         lenseq = ceil(log10(idxtile));
%         st = [0 0 0];
%         bin = chunk_size_ijk'*2.^(lenseq-1:-1:0);
%         indxyz = zeros(lenseq,3);
%         for iseq = 1:lenseq
%             is = seq(iseq);
%             temp = fliplr(dec2bin(str2num(is)-1,3)); %#ok<ST2NM>  % coordinates wrto xyz
%             for ii=1:3 % x y z
%                 indxyz(iseq,ii) = str2num(temp(ii));  %#ok<ST2NM>
%                 st(ii) = st(ii) + bin(ii,iseq)*indxyz(iseq,ii);
%             end
%         end
%         RR(laden_chunk_index,:)=[st st+chunk_size_ijk];
%         
%         octree_path = arrayfun(@str2double, seq) ;
%         chunk_ijk1 = chunk_ijk1s_from_octree_paths(octree_path) ;
%         st_check = (chunk_ijk1-1) .* chunk_size_ijk ;
%         assert(isequal(st, st_check)) ; 
%     end

    octree_paths = zeros(laden_chunk_count, zoom_step_count) ;
    parfor laden_chunk_index = 1:laden_chunk_count ,
        laden_chunk_file_path = laden_chunk_file_paths{laden_chunk_index} ;
        path_elements = strsplit(laden_chunk_file_path, filesep) ;
        octree_path_as_string = [path_elements{end-zoom_step_count:end-1}] ;  % octree path as an array of chars, each char a digit on 1...8
        octree_path = arrayfun(@str2double, octree_path_as_string) ;
        octree_paths(laden_chunk_index,:) = octree_path ;
    end
    chunk_ijk1s = chunk_ijk1s_from_octree_paths(octree_paths) ;
    region_lower_corners_ijk0 = (chunk_ijk1s-1) .* chunk_size_ijk ;
    RR = [region_lower_corners_ijk0 region_lower_corners_ijk0+chunk_size_ijk] ;
    
    %assert(isequal(RR, RR_check)) ;    
end

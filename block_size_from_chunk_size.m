function block_size_ijk = block_size_from_chunk_size(chunk_size_ijk)
    % This divides up the chunk_size in each dim as needed to get all
    % chunks <= 128 voxels on a size.
    block_size_ijk = chunk_size_ijk/2 ;
    max_block_size_ijk = [128 128 128] ;
    while any( block_size_ijk > max_block_size_ijk ) ,
        is_dim_too_big = (block_size_ijk > max_block_size_ijk) ;
        dim_divisor = is_dim_too_big + 1 ;  % 2 if too big, 1 if OK
        block_size_ijk = block_size_ijk ./ dim_divisor ;
    end
end

function mergeh5_chunks(sample_date, mask_threshold)
    input_folder_path = sprintf('/nrs/mouselight/SAMPLES/%s-prob', sample_date) ;
    output_folder_path = sprintf('/groups/mousebrainmicro/mousebrainmicro/cluster/Reconstructions/%s/whole-brain-p-map-as-h5', sample_date) ;
    mergeh5_chunks_given_paths(output_folder_path, input_folder_path, mask_threshold) ;
end



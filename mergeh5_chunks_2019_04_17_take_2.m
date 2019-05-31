this_file_path = mfilename('fullpath') ;
this_folder_path = fileparts(this_file_path) ;
configuration_file_path = fullfile(this_folder_path, 'config_files', '2019-04-17-mergeh5-take-2.cfg') ;
mergeh5_chunks(configuration_file_path) ;

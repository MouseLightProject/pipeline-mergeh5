function modpath()
    this_script_path = mfilename('fullpath') ;
    this_script_folder_path = fileparts(this_script_path) ;        
    addpath(genpath(fullfile(this_script_folder_path, 'common'))) ;
    addpath(genpath(fullfile(this_script_folder_path, 'mouselight_toolbox'))) ;
end

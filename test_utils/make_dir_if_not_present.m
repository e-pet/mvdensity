function make_dir_if_not_present(filename)

[pathstr, ~, ~] = fileparts(filename);
    if ~exist(pathstr,'dir') 
        mkdir(pathstr);
    end
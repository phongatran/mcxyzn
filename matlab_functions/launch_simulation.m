function [] = launch_simulation(cfg)
    if (ispc)
        bin_name = 'mcxyzn.exe';
    elseif (ismac)
        bin_name = 'mcxyzn.mac';
    elseif (isunix)
        bin_name = 'mcxyzn.linux';
    else
        fprintf('Could not find the appropriate binary \n');
    end
    filepath = which(bin_name);
    system_command_string = [sprintf('%s',filepath),' ',cfg.name];

    [status] = system(system_command_string);
end

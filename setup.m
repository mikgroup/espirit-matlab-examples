% To begin, add '$(TOOLBOX_PATH)/matlab' to the library path. The
% environment variable TOOLBOX_PATH needs to be set to the base
% directory of the reconstruction tools package.

if isempty(getenv('TOOLBOX_PATH'))
	error('Environment variable TOOLBOX_PATH is not set.');
end

addpath(strcat(getenv('TOOLBOX_PATH'), '/matlab'));


% print version information
bart('version -V');


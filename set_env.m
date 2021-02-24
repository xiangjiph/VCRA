% Setup environments
% clc;close all; clear; 
[~, computer_name] = system('hostname');
switch computer_name(1:end-1)
    case ''        
        SCRIPT_PATH = '';
    otherwise
        SCRIPT_PATH = '.';
        disp('Computer not recognized. Use the current directory.');
end
addpath(genpath(SCRIPT_PATH));
clear;clc;
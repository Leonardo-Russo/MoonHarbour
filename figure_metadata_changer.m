%% Figure Metadata Changer

close all
clear
clc

%% Core

addpath("Output/Stuff/");

filepath = "reference_trajectory_top";

% Load the .fig file
fig = openfig(strcat("Output/Stuff/", filepath, ".fig"));



% Save the figure as a PNG with 600 DPI
savefig(fig, strcat("Output/Stuff/", filepath, ".fig"));
print(fig, strcat("Output/Stuff/", filepath, ".png"), '-dpng', '-r600');          % 300 DPI

% Close the figure
close(fig);

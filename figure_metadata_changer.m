%% Figure Metadata Changer

close all
clear
clc

%% Core

addpath("Stuff/");

filepath = "saturated_control_90m";

% Load the .fig file
fig = openfig(strcat("Stuff/", filepath, ".fig"));

title('');

% Save the figure as a PNG with 600 DPI
savefig(fig, strcat("Stuff/", filepath, ".fig"));
print(fig, strcat("Stuff/", filepath, ".png"), '-dpng', '-r600');          % 300 DPI

% Close the figure
close(fig);

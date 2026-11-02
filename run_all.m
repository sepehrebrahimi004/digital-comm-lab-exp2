%RUN_ALL  Digital Communications Lab - Experiment 2: source and channel coding.
%   Writes every figure and the recovered texts to ./results.
%   Requires the Communications Toolbox and the Image Processing Toolbox.

clear; close all; clc;
addpath(fullfile(fileparts(mfilename('fullpath')), 'src'));
rng(1405);

q2_1_huffman_qam;   clearvars
q2_2_ldpc_image;    clearvars
q2_3_huffman_rs;    clearvars
close all;
disp('All results written to ./results');

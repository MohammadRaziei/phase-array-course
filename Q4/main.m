%% Q4
clear; close all; clc;
mkdir results
addpath ../common/

sll = 20;
run_for_sll(sll);

sll = 30;
run_for_sll(sll);

sll = 40;
run_for_sll(sll);

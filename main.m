
clear all;
clc;
close all;
mex cec13_func.cpp -DWINDOWS
%  inputs:
% N:  Number of agents.
% max_it: Maximum number of iterations (T).
% ElitistCheck: If ElitistCheck=1, algorithm runs with modified Kbest and
% if =0, runs with Kbest = N.
% Rpower: power of 'R'
% F_index: The index of the test function. See tables I of the mentioned article.
%          Insert your own objective function with a new F_index in 'test_functions.m'
%          and 'test_functions_range.m'.
%            For demonstration, Rosenbrock’s Function is taken as example

%  outputs:22
% Fbest: Best result. 
% Lbest: Best solution. The location of Fbest in search space.
% BestChart: The best so far Chart over iterations. 
% MeanChart: The average fitnesses Chart over iterations.
global C_dia min_pres P_length nnode npipe y;

for k=1:1
    
inpname='balerma.inp';     % Hanoi Network Benchmark Problem

epanet(inpname);

%%
% System information
njunc=getdata('EN_NODECOUNT');    % Total junction number
ntank=getdata('EN_TANKCOUNT');    % Reservoir/tank number
nnode=njunc-ntank;                % Node number
npipe=getdata('EN_LINKCOUNT');    % Pipe number
P_length=getdata('EN_LENGTH')';
y='x';
nvars=npipe;
%%
% commercial pipe diameter
C_dia=[304.8 406.4 508 609.6 762 1016];
N_dia=length(C_dia);
min_pres=30;
 N=50; 
 max_it=600; 
 ElitistCheck=1; Rpower=1;
 min_flag=1; % 1: minimization, 0: maximization

 F_index=10; 

 %%
[Fbest,Lbest,BestChart,MeanChart]=CLGSA(F_index,N,max_it,ElitistCheck,min_flag,Rpower)
 %%


min_cost=min(Fbest);
max_cost=max(Fbest);
average_cost=mean(Fbest);
SD_cost=std(Fbest);
end

 figure(1);
 plot(BestChart,'--k');
 
title(['\fontsize{12}\bf F',num2str(F_index)]);
xlabel('\fontsize{12}\bf Iteration');ylabel('\fontsize{12}\bf Best-so-far');
legend('\fontsize{10}\bf CLGSA','NorthEast');
disp(['k= ',num2str(k),'   FF= ',num2str(Fbest)])
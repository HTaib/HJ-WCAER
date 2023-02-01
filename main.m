%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is part of research works on Data clustering using HJ-WCAER                             %
% Paper co-author and programmer: Hasnanizan Taib                                                    %
%                                                                                                                                    %
% Email: hasnanizan.taib@gmail.com                                                                             %
%           hasnanizan.taib@utb.edu.bn                                                                             %
%                                                                                                                                    %
% Main paper to be cited:                                                                                                %
% Data clustering using hybrid water cycle algorithm and a local pattern search method %
% by H. Taib & A. Bahreininejad                                                                                      %
% Advances in Engineering Software                                                                               %
% DOI: https://doi.org/10.1016/j.advengsoft.2020.102961                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;clc;close all;

format long g
format compact
rng('default')
for k=1:1 % Number of independent runs
    %%%Choose the benchmark dataset (Uncomment one of them)
    %%%The Iris dataset--------------------------------------------------------
    data=load('iris');        
    X=data.iris;              
    nClust = 3;     % initial number of cluster centres, for benchmark dataset, they were known
    %%%The Glass dataset--------------------------------------------------------
    % data=load('Glass');        
    % X=data.Glass;              
    % nClust = 6;     % initial number of cluster centres, for benchmark dataset, they were known
    %%%The Wisconsin Breast Cancer dataset--------------------------------------------------------
    % data=load('wbc');        
    % X=data.wbc2;              
    % nClust = 2;     % initial number of cluster centres, for benchmark dataset, they were known
    %%%--------------------------------------------------------
    nVars=size(X,2);      % Number of features or variables
    %%%Choose the objective function (Uncomment one of them)
    objective_function=@(s) Euclidean(s, X);     % Euclidean distance
    % objective_function=@(s) DBIndex(s, X);        %DB index
    %%%--------------------------------------------------------
    constraints=@const;
    LB= repmat(min(X),nClust,1);   % Lower Bound of each column of variables
    UB= repmat(max(X),nClust,1);  % Upper Bound of each column of variables
    %%%Experimental setting
    max_it=1000;      % Maximum number of generations/iterations
    nPop=50;           % Size of population 
    
    %%%Algorithm parameter setting
    Nsr=12;               %Number of rivers + sea
    dmax=1e-8;         %Evaporation condition constant
    
    [Xmin,Fmin,SUM_Constraints,NFEs,Elapsed_Time,iters]=HJWCAER(objective_function,constraints,LB,UB,nClust,nVars,nPop,Nsr,dmax,max_it)

    disp(['Run: ',num2str(k),'   Fmin= ',num2str(Fmin),'  Summation Constraint Violations:  ',num2str(SUM_Constraints)]);
    fprintf('\nX values ='); fprintf('\t%12.15f   ', Xmin); fprintf('\n');
    F(k)=Fmin;
end

[MinimumFunctionValues index]=min(F)
AverageFunctionValues=mean(F)
MaximumFunctionValues=max(F)
StandardDeviation=std(F)
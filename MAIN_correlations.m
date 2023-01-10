% Demonstrate algorithm and results from an alternate fragility correlations approach
% Created by Jack Baker, 4/20/2021, last updated 1/10/2023
%
% The algorithm demonstrated here is described in the following document,
% and the results from this calculation are the "Simple Example" presented
% in this manuscript.
%
% Jack W. Baker, Ed Almeter, Dustin Cook, Abbie Liel, and Curt Haselton
% (2023) "A model for partially dependent component damage fragilities in 
% seismic risk analysis," (in review).


clear; close all; clc;


%% Specify input parameters

nSims = 10000; % how many Monte Carlo simulations to run
makeFigs = 1; % = 1 to plot results figures

% EDPs (generic EDPs--for example, SDRs on nEDP floors)
edp.nEDP = 4; % how many EDPs
edp.beta = 0.5; % log standard deviations of EDPs
edp.rhoEDP = 0.6; % correlation between EDPs

% components (assume a single damage state, for simplicity)
comp.nCompTypes = 5;
edp.theta = 0.02; % median EDP value
comp.theta = 0.05;
comp.beta = 0.5;
comp.quantity = 1; % number of components per EDP

%% build component inventory -- assume every component has the same characteristics for simplicity, but this could be modified
for i = 1:comp.nCompTypes
    component(i).theta = comp.theta;
    component(i).beta = comp.beta;
    component(i).quantity = comp.quantity*ones(1,edp.nEDP); % entry j is quantity associated with EDP j
end

%% some simple calculations to check the reasonableness of the above parameters
z = (log(edp.theta) - log(comp.theta)) / sqrt(edp.beta^2 + comp.beta^2) % how far is the the median EDP value from the median component capacity 
pf = normcdf( z) % what is the component probability of failure
E_N = pf*edp.nEDP*comp.nCompTypes % what is the expected number of damaged components

%% Run the Monte Carlo analyses
results = fn_simulate_damage(nSims, edp, component);



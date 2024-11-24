clear; clc; close all;
%test the program by a SDoF problem
%by J.S. Yang
%date: 2019-11-05
%representative point set
%date: 2021-06-23
%update the call of PDEM_solve function
omega = (5*pi/4:pi/200:7*pi/4)';
asgn_prob = [1/200, 1/100*ones(1, 99), 1/200]';
np = length(omega);
%system response
x0 = 0.1;
v0 = 0;
dt = 0.001;
t0 = (0:dt:10)';
nt = length(t0);
d=zeros(length(t0),np);
v=zeros(length(t0),np);
for ii = 1:1:np
    d(:,ii) =  x0*cos(omega(ii)*t0);
    v(:,ii) = -x0*omega(ii)*sin(omega(ii)*t0);
end
thres = inf;
dt_ratio = 1;
[tm, rm, prob] = PDEM_solve(asgn_prob, d, v, t0, dt, dt_ratio, thres);
%%
PDEM_post(tm, rm, prob, 'surf', [0.9,1.1], [-0.2,0.2], 'Displacement[m]');
PDEM_post(tm, rm, prob, 'contour', [0.9,1.1], [-0.2,0.2], 'Displacement[m]');
PDEM_post(tm, rm, prob, 'curve', [0.9,1.0,1.1], [-0.2,0.2], 'Displacement[m]');

%% References:
% -------------------------------------------------------------------------------
% Li J, Chen JB. (2009) Stochastic dynamics of structures. Wiley, Singapore.
% Li J, Chen JB. (2004) Probability density evolution method for dynamic response analysis of structures with uncertain parameters. Comput Mech 34(5): 400-409. http://link.springer.com/10.1007/s00466-004-0583-8
% Chen JB, Yang JS, Jensen HA. (2020) Structural optimization considering dynamic reliability constraints via probability density evolution method and change of probability measure. Struct Multidisc Optim 62: 2499–2516. https://doi.org/10.1007/s00158-020-02621-4
% Yang JS, Jensen HA, Chen JB. (2022) Structural optimization under dynamic reliability constraints utilizing probability density evolution method and metamodels in augmented input space. Struct Multidisc Optim 65: 107. https://doi.org/10.1007/s00158-022-03188-y
% Yang JS, Chen JB, Jensen HA. (2022) Structural design optimization under dynamic reliability constraints based on the probability density evolution method and highly-efficient sensitivity analysis. Probab Eng Mech 68: 103205. https://doi.org/10.1016/j.probengmech.2022.103205

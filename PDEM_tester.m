clear; clc; close all;
%test the program by a SDoF problem
%by J.S. Yang
%date: 2019-11-05
%representative point set
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
[tm, rm, prob] = PDEM_solve(asgn_prob, d, v, t0, dt, thres);
%%
PDEM_post(tm, rm, prob, 'surf', [0.9,1.1], [-0.2,0.2], 'Displacement[m]');
PDEM_post(tm, rm, prob, 'contour', [0.9,1.1], [-0.2,0.2], 'Displacement[m]');
PDEM_post(tm, rm, prob, 'curve', [0.9,1.0,1.1], [-0.2,0.2], 'Displacement[m]');
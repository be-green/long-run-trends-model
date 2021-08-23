function [phi,alpha,lambda,mu,hc_loss,n_periods,g,delta,omega,sigma,rho,v,p_z,kappa,theta_grid,theta0,xi_constant] = ...
          parse_model_params_v2(paramvec,H_inside,n_gridpoints)
% [phi,alpha,lambda,mu,hc_loss,n_periods,g,delta,omega,sigma,rho,v,p_z,kappa,theta_grid,theta0] = ...
%  parse_model_params_v1(paramvec,H_inside,n_gridpoints)
% INPUTS:
%   - paramvec: list of parameters, which includes
%        (1): phi, the rate at which skilled workers move up the grid
%        (2): alpha, the probability of being displaced, given an
%             innovation arrives
%        (3): expected_hc_loss, d * alpha x omega (arrival rate of new innovations)
%        (4): CES parameter on L - CES parameter on H
%        (5): CES parameter on L 
%        (6): kappa, size of jump in xi
%        (7): a, log odds ratio capturing higher probability of exposed
%             workers being displaced
%        (8): mu, parameter on labor in outer nest
%        (9): lambda, parameter on labor in inner nest
%       (10): xi_constant
%   - H_inside: dummy to indicate whether H is in the inner nest
%   - n_gridpoints: number of gridpoints to use for theta grid 
%
% OUTPUTS:
%   - phi, the rate at which skilled workers move up the grid
%   - alpha, the probability of being displaced, given an
%            innovation arrives
%   - lambda, parameter on labor in inner nest
%   - mu, parameter on labor in outer nest
%   - hc_loss, workers lose hc_loss*theta if displaced
%   - n_periods, number of periods in a 5 year interval (fixed inside script)
%   - g, population growth (fixed inside script)
%   - delta, death rate (fixed inside script)
%   - omega, arrival rate of new innovations
%   - sigma, outer nest CES parameter
%   - rho, inner nest CES parameter
%   - v, DRS parameter (Fixed at 1 currently)
%   - p_z, probability of displacement for exposed workers post shock
%   - kappa, size of shock to xi
%   - theta_grid, support for human capital levels, on the unit interval
%   - theta0: bottom of theta grid
%   - xi_constant: intercept term to add for drift of xi in VAR
      
paramvec = paramvec';

% Ex-ante fixed parameters
% 0 is geometric increasing from theta
% 1 is geometric decreasing from 1
theta_order = 0;

% unconditional growth

% number of iterations associated with a single shock
% since a shock is 5 years, this corresponds to there being
% 4 iterations of the VAR a year
n_periods = 60;

g = exp(log(0.02 + 1) / (n_periods / 5)) - 1;

% death rate
% corresponds to an average 40 year working life
% also, newborn agents start at the bottom of the grid, so it's like
% effective death rate is g% higher;
delta =  exp(log(0.025 + 0.02 + 1) / (n_periods / 5)) - 1;

% set up variables in accordance w/ our model

% rate of moving up the ladder
phi = (paramvec(1,:));
% phi = paramvec(1,:);

% conditional fall probability
alpha = (paramvec(2,:));


% arrival rate of technology
% in this case this is also the probability
% at period t
% of state 1 happening in period t + 1
% these shocks are IID so this is true for both
% initial states
omega = 5/n_periods; %fixing this ex ante in this round (1x / year)

d_x_omega_x_alpha = (paramvec(3,:));
d = d_x_omega_x_alpha / (alpha * omega);
hc_loss = 1-exp(-d);
% omega = paramvec(3,:);

% sigma is outer nest exponent
% rho is inner nest exponent
if H_inside == 1
    sigma = (paramvec(5,:));
    rho = sigma - (paramvec(4,:));
else
    rho = (paramvec(5,:));
    sigma = rho - (paramvec(4,:));
end

% loss in human capital. Right now this is exp(-d) times prior theta level
kappa = paramvec(6,:);

% should this be paramvec(7:,)? ;
p_z = exp(paramvec(7,:) + log(alpha)) / (1 + exp(paramvec(7,:) + log(alpha)));

% CES share parameters
lambda = paramvec(9,:); % inner nest
mu = paramvec(8,:); % outer nest

% DRS parameter. Fixed at the start
v = 1;

% intercept on xi in the VAR;
xi_constant = paramvec(10,:);

% setting up theta grid
theta0 = 0.05;

if theta_order == 0
    % growth_rate = paramvec(8,:);
    % n_gridpoints = floor(-log(theta0) / log(1 + growth_rate));
    growth_rate = exp((-log(theta0)) / n_gridpoints) - 1;
    theta_grid = (theta0).*((1 + growth_rate).^(1:(n_gridpoints))); 
else 
    growth_rate = exp((log(theta0)) / n_gridpoints) - 1;
    theta_grid = 1 - (1.*((1 + growth_rate).^(1:(n_gridpoints))));
end

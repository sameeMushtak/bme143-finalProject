% BME 143 - Final Model Implementation
% pdeScript generates N simulated Drosophila embryos and stores the
% simulations in a .mat file
%% Setup
clear
close all
% Symmetry Parameter, 1-D Cartesian
m = 0;
% Variables determining what simulation to run
J_B_fac = 1; % 0.5,1,1.5,2
nu_fac = 1; % 1e-3,1e-2,0.1,0.2,0.5,1
eta_fac = 1; % 0,1
staufen_frac = 1; % Wild-Type: 1, Mutant: 0.93
x_center = 17.5; % Default 17.5
x_width = 12.5; % Default 12.5
delta_radius = 5;
filename = 'default.mat';
% Number of simulations
N = 100;
% Number of position points
x_pts = 500;
% Number of time points
t_pts = 500;
% Time Vector, Steady State assumed to occur after one hour (see page 2)
t = linspace(0,4000,t_pts);
% Preallocation
x  = zeros(x_pts,N);
EL = zeros(N,1);
mu = zeros(N,1);
M  = zeros(x_pts,N,6);

%% Solution of System of PDEs
for i = 1:N
    % Position Vector, Length of embryo in um
    EL_mean = 460;
    % 10% fluctuations contained within 3 standard deviations
    EL(i) = (EL_mean/30)*randn(1)+EL_mean;
    x(:,i) = linspace(0,EL(i),x_pts);
    % Random parameters sampled from Gaussian distribution
    J_B = p_gauss(0.08*J_B_fac,0.05);
    J_R = p_gauss(0.08,0.05);
    mu_hb = p_gauss(4e-3,0.3);
    alpha = p_gauss(6e-3,0.3);
    beta = p_gauss(4e-3,0.3);
    K = p_gauss(4.5e-2,0.3);
    % Implementation of correlation between J_B and J_R
%     mean_J = [0.08*J_B_fac 0.08];
%     rho_J = 1;
%     p_J = 0.45;
%     Sigma_J = [(mean_J(1)*p_J)^2 rho_J*(mean_J(1)*p_J)*(mean_J(2)*p_J); rho_J*(mean_J(1)*p_J)*(mean_J(2)*p_J) (mean_J(2)*p_J)^2];
%     pair_J = mvnrnd(mean_J,Sigma_J);
%     while (pair_J(1) < 0) || (pair_J(2) < 0)
%         pair_J = mvnrnd(mean_J,Sigma_J);
%     end
%     J_B = pair_J(1);
%     J_R = pair_J(2);
    % Generating correlated mu and mu_Hb values
    mean_mu = [7e-4 6e-3];
    rho_mu = 1;
    p_mu = 0.3;
    Sigma_mu = [(mean_mu(1)*p_mu)^2 rho_mu*(mean_mu(1)*p_mu)*(mean_mu(2)*p_mu); rho_mu*(mean_mu(1)*p_mu)*(mean_mu(2)*p_mu) (mean_mu(2)*p_mu)^2];
    pair_mu = mvnrnd(mean_mu,Sigma_mu);
    while (pair_mu(1) < 0) || (pair_mu(2) < 0)
        % This messes up the correlation, but it should not have a major
        % impact
        pair_mu = mvnrnd(mean_mu,Sigma_mu);
    end
    mu(i) = pair_mu(1);
    mu_Hb = pair_mu(2);
    
    % Location of synthesis of Bcd and repressor protein
    x_B = (x_center-x_width)+x_width*rand(1);
    x_R = EL(i)-((x_center-x_width)+x_width*rand(1));
    
    nu = nu_fac*3;
    eta = eta_fac*1.0;
    D = 10;
    D_hb = 0.5;
    D_Hb = 1;
    gamma = 0.01;
    
    % frac is a random number between staufen_frac and 1
    frac = staufen_frac + (1-staufen_frac)*rand(1);
    params = [J_B J_R mu(i) mu_hb mu_Hb alpha beta K x_B x_R nu eta...
        D D_hb D_Hb gamma frac delta_radius];
    sol = pdepe(m,@(x,t,u,dudx) pdefun(x,t,u,dudx,params),@pdeic,@pdebc,x(:,i),t);
    % Add a plane for [B]+[C]
    sol(:,:,6) = sol(:,:,1)+sol(:,:,3);
    M(:,i,:) = sol(end,:,:);
    X = sprintf('%i/%i',i,N);
    disp(X)
end

save(filename,'N','x','EL','mu','M')

% Generates a random number from a Gaussian distribution with mean mu and
% standard deviation p_sigma*mu
% Rerolls if the generated value is negative
function obs = p_gauss(mu,p_sigma)
    obs = -1;
    while obs < 0
        obs = (p_sigma*mu)*randn(1)+mu;
    end
end
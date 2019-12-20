% BME 143 - Model Extension
%% Setup
clear
close all
% Symmetry Parameter, 1-D Cartesian
m = 0;
% Diffusion coefficient for Bcd is fixed at 10
D_B = 10;
% Diffusion coefficients for corepressor and Bcd-corepressor complex are
% varied with respect to Bcd's
D_R = 10;
D_C = 10;
filename = 'simulations_(1,1).mat';
% Number of simulations
N = 100;
% Number of position points
x_pts = 250;
% Number of time points
t_pts = 250;
% Time Vector, Steady State assumed to occur after one hour
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
    % Random parameters with Gaussian distribution
    J_B = p_gauss(0.08,0.05);
    J_R = p_gauss(0.08,0.05);

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
    
    mu_hb = p_gauss(4e-3,0.3);
    alpha = p_gauss(6e-3,0.3);
    beta = p_gauss(4e-3,0.3);
    K = p_gauss(4.5e-2,0.3);
    
    x_B = 5+25*rand(1);
    x_R = EL(i)-(5+25*rand(1));
    nu = 3;
    eta = 1.0;
    D_hb = 0.5;
    D_Hb = 1;
    gamma = 0.01;
    
    params = [J_B J_R mu(i) mu_hb mu_Hb alpha beta K x_B x_R nu eta...
        D_B D_R D_C D_hb D_Hb gamma];
    sol = pdepe(m,@(x,t,u,dudx) pdefunExtended(x,t,u,dudx,params),@pdeic,@pdebc,x(:,i),t);
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
% Pre-Figure 2 Calculations
%% Setup - Single Test
clear
close all
test = load('fig2c_mats/simulations_staufen_mutant.mat');

%% Setup - J_B multiplier tests
clear
close all
test(1) = load('simulations_Bcd_dose_1-2.mat');
test(2) = load('figS1_mats/simulations_delta-5.mat');
test(3) = load('simulations_Bcd_dose_3-2.mat');
test(4) = load('simulations_Bcd_dose_2-1.mat');

%% Setup - Bcd/corepressor (nu) binding tests
clear
close all
test(1) = load('fig2a_mats/simulations_Bcd(A52-56)_mutant.mat');
test(2) = load('fig2a_mats/simulations_nu_0.03.mat');
test(3) = load('fig2a_mats/simulations_nu_0.3.mat');
test(4) = load('fig2a_mats/simulations_nu_0.6.mat');
test(5) = load('fig2a_mats/simulations_nu_1.5.mat');
test(6) = load('figS1_mats/simulations_delta-5.mat');

%% Setup - Delta Radius tests
clear
close all
for i = 1:30
    filename = sprintf('figS1_mats/simulations_delta-%i.mat',i);
    test(i) = load(filename);
end

%% Setup - x_width tests
clear
close all
test(1) = load('figS2_mats/simulations_17.5pm00.0.mat');
test(2) = load('figS2_mats/simulations_17.5pm02.5.mat');
test(3) = load('figS2_mats/simulations_17.5pm05.0.mat');
test(4) = load('figS2_mats/simulations_17.5pm07.5.mat');
test(5) = load('figS2_mats/simulations_17.5pm10.0.mat');
test(6) = load('figS1_mats/simulations_delta-5.mat');
test(7) = load('figS2_mats/simulations_17.5pm15.0.mat');
test(8) = load('figS2_mats/simulations_17.5pm17.5.mat');

%% Setup - x_center tests
clear
close all
test(1) = load('figS3_mats/simulations_12.5pm12.5.mat');
test(2) = load('figS3_mats/simulations_15.0pm12.5.mat');
test(3) = load('figS1_mats/simulations_delta-5.mat');
test(4) = load('figS3_mats/simulations_20.0pm12.5.mat');
test(5) = load('figS3_mats/simulations_22.5pm12.5.mat');

%% Calcualations
x_Hb = zeros(test(1).N,length(test));
for i = 1:length(test)
    x_Hb(:,i) = x_boundary(test(i),5,0.5)./test(i).EL;
end
% Compute mean and standard deviation in position of Hb boundary
stats_x_Hb = [mean(x_Hb); std(x_Hb)];
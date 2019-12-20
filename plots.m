% Plots
%% Setup
clear
close all
fig1 = load('figS1_mats/simulations_delta-5.mat');
fig2a(1) = load('fig2a_mats/simulations_Bcd(A52-56)_mutant.mat');
fig2a(2) = load('fig2a_mats/simulations_nu_0.03.mat');
fig2a(3) = load('fig2a_mats/simulations_nu_0.3.mat');
fig2a(4) = load('fig2a_mats/simulations_nu_0.6.mat');
fig2a(5) = load('fig2a_mats/simulations_nu_1.5.mat');
fig2a(6) = load('figS1_mats/simulations_delta-5.mat');
fig2b = load('fig2b_mats/simulations_defective_Hb.mat');
fig2c = load('fig2c_mats/simulations_staufen_mutant.mat');
for i = 1:30
    filename = sprintf('figS1_mats/simulations_delta-%i.mat',i);
    figS1(i) = load(filename);
end
figS2(1) = load('figS2_mats/simulations_17.5pm00.0.mat');
figS2(2) = load('figS2_mats/simulations_17.5pm02.5.mat');
figS2(3) = load('figS2_mats/simulations_17.5pm05.0.mat');
figS2(4) = load('figS2_mats/simulations_17.5pm07.5.mat');
figS2(5) = load('figS2_mats/simulations_17.5pm10.0.mat');
figS2(6) = load('figS1_mats/simulations_delta-5.mat');
figS2(7) = load('figS2_mats/simulations_17.5pm15.0.mat');
figS2(8) = load('figS2_mats/simulations_17.5pm17.5.mat');
figS3(1) = load('figS3_mats/simulations_12.5pm12.5.mat');
figS3(2) = load('figS3_mats/simulations_15.0pm12.5.mat');
figS3(3) = load('figS1_mats/simulations_delta-5.mat');
figS3(4) = load('figS3_mats/simulations_20.0pm12.5.mat');
figS3(5) = load('figS3_mats/simulations_22.5pm12.5.mat');
figS4 = load('extensionStats.mat');
%% Figure 1b
figure
hold on
profile_plot(fig1,6,'b')
xlabel('{\itx}/EL')
ylabel('Normalized Density')
title('Mean Density Profile of Total Bcd Protein')
saveas(gcf,'pngs/fig1b.png')

%% Figure 1c
figure
hold on
profile_plot(fig1,1,'b')
profile_plot(fig1,2,'r')
xlabel('{\itx}/EL')
ylabel('Normalized Density')
text(0.14,0.8,'Unbound Bcd')
text(0.57,0.8,'Unbound Corepressor')
title('Mean Density Profile of Unbound Bcd and Unbound Corepressor Protein')
saveas(gcf,'pngs/fig1c.png')

%% Figure 1d
figure
hold on
profile_plot(fig1,5,'b')
xlabel('{\itx}/EL')
ylabel('Normalized Density')
title('Mean Density Profile of Hb Protein')
text(0.1,0.5,'{\itx}_{Hb} = 0.49(\pm0.02)EL')
saveas(gcf,'pngs/fig1d.png')

%% Figure 1e
lambda = sqrt(10 ./ fig1.mu);
figure
histogram(lambda)
xlabel('\lambda (\mum)')
ylabel('Number of Embryos')
title('Histogram of length scale describing Total Bcd density profile')
saveas(gcf,'pngs/fig1e.png')

%% Figure 1f
% Find position closest to normalized density of 0.17 for Bcd, 0.5 for Hb
x_Bcd = x_boundary(fig1,6,0.17);
x_Hb = x_boundary(fig1,5,0.5);
figure
hold on
histogram(x_Bcd./fig1.EL)
histogram(x_Hb./fig1.EL)
xlabel('{\itx}_{Bcd} / EL; {\itx}_{Hb} / EL')
ylabel('Number of Embryos')
legend('Bcd','Hb')
title('Histograms of Positions of Bcd and Hb Boundaries')
saveas(gcf,'pngs/fig1f.png')

%% Figure 1g
figure
hold on
plot(fig1.EL,x_Bcd,'ko')
plot(fig1.EL,x_Hb,'gs')
fplot(@(t) (fig1.EL \ x_Hb)*t,'b')
xlabel('EL (\mum)')
ylabel('Boundary Position (\mum)')
legend('{\itx}_{Bcd}','{\itx}_{Hb}')
xlim([420 500])
ylim([150 300])
title('Positions of Bcd and Hb Boundaries vs. Embryo Length')
text(487,235,'{\itx}_{Hb}=0.49EL')
saveas(gcf,'pngs/fig1g.png')

%% Figure 2a - Main
figure
hold on
profile_plot(fig2a(1),4,'b')
xlabel('{\itx}/EL')
ylabel('Normalized Density')
text(0.1,0.5,'{\itx}_{\ithb} = 0.77(\pm0.10)EL')
title('Mean Density Profile of hb mRNA in Bcd(A52-56) mutants')
saveas(gcf,'pngs/fig2a-Main.png')

%% Figure 2a - Inset
x_hb = zeros(fig2a(1).N,length(fig2a));
for i = 1:length(fig2a)
    x_hb(:,i) = x_boundary(fig2a(i),4,0.5)./fig2a(i).EL;
end
mean_x_hb = mean(x_hb);
std_x_hb = std(x_hb);
figure
hold on
errorbar([0.003 0.03 0.3 0.6 1.5 3],mean_x_hb,std_x_hb,'ko-')
xlabel('\nu (\mum\cdots^{-1})')
ylabel('{\itx}_{hb} / EL')
set(gca, 'XScale', 'log')
xlim([1e-3 5])
ylim([0.4 1])
title('Position of {\ithb} boundary in Bcd(A52-56)-like mutants')
saveas(gcf,'pngs/fig2a-Insert.png')

%% Figure 2b
figure
hold on
profile_plot(fig2b,5,'b')
xlabel('{\itx}/EL')
ylabel('Normalized Density')
text(0.1,0.5,'{\itx}_{Hb} = 0.46(\pm0.02)EL')
title('Mean Density Profile of Hb Protein for \eta=0')
saveas(gcf,'pngs/fig2b.png')

%% Figure 2c
figure
hold on
profile_plot(fig2c,5,'b')
xlabel('{\itx}/EL')
ylabel('Normalized Density')
text(0.04,0.5,'{\itx}_{Hb} = 0.43(\pm0.04)EL')
title('Mean Density Profile of Hb Protein in {\itstaufen} mutants')
saveas(gcf,'pngs/fig2c.png')

%% Figure 3
figure
hold on
x1 = 0:0.15:0.45;
x2 = 0:0.05:0.45;
plot(x1,[0.0159 0.0205 0.0219 0.0268],'m>-')
plot(x1,[0.0219 0.0209 0.0219 0.0277],'cv-')
plot(x1,[0.0207 0.0194 0.0219 0.0294],'y<-')
plot(x1,[0.0223 0.0192 0.0219 0.0243],'b^-')
plot(x1,[0.0203 0.0184 0.0219 0.0250],'gd-')
plot(x1,[0.0195 0.0255 0.0270 0.0326],'rs-')
plot(x2,[0.0203 0.0219 0.0274 0.0325 0.0348 0.0520 0.0604 0.0741 0.0935 0.1014],'ko-')
title('Variation in Hb boundary vs. Variation in parameters')
xlabel('relative magnitude of fluctuations')
ylabel('\sigma_{{\itx}_{Hb}}')
saveas(gcf,'pngs/fig3.png')

%% Figure S1 - Varying Delta_Radius
figure
hold on
stat_plot(figS1,1:length(figS1))
xlabel('Radius of Dirac delta (\mum)')
ylabel('{\itx}_{Hb}/EL')
title('Position of Hb boundary as a function of Dirac delta radius')
saveas(gcf,'pngs/figS1.png')

%% Figure S2 - Varying x_width
figure
hold on
stat_plot(figS2,0:5:35)
xlabel('Diameter of possible x locations (\mum)')
ylabel('{\itx}_{Hb}/EL')
title('Position of Hb boundary as a function of range of source locations')
saveas(gcf,'pngs/figS2.png')

%% Figure S3 - Varying x_center
figure
hold on
stat_plot(figS3,12.5:2.5:22.5)
xlabel('Diameter of possible x locations (\mum)')
ylabel('{\itx}_{Hb}/EL')
title('Position of Hb boundary as a function of range of source locations')
saveas(gcf,'pngs/figS3.png')

%% Figure S4 - Altering Diffusion Coefficients
[X,Y] = meshgrid(figS4.D,figS4.D);
figure
hold on
surf(X,Y,figS4.mean_x_Hb)
xlabel('{\itD}_C/{\itD}_B')
ylabel('{\itD}_R/{\itD}_B')
zlabel('{\itx}_{Hb}')
title('Position of Hb boundary vs. Diffusion Coefficient')

%% Plotting Functions
function profile_plot(s,cmpd,plot_options)
% Plots 10 sample curves with error bars
% Inputs
    % s - Struct generated by loading .mat file created by pdeScript.m
    % or modelExtension.m
    % cmpd - Compound identifier. Integer from 1 to 6.
        % cmpd = 1 ~ [B]
        % cmpd = 2 ~ [R]
        % cmpd = 3 ~ [C]
        % cmpd = 4 ~ [hb]
        % cmpd = 5 ~ [Hb]
        % cmpd = 6 ~ [B]+[C]
    % plot_options - String containing information for plot modifications
    
    rho = normalize_density(s,cmpd);
    % for i = 1:s.N
    % Plot 10 normalized density profiles
    for i = 1:10
        plot(s.x(:,i)/s.EL(i),rho(:,i),plot_options)
    end
    x_rho = zeros(s.N,5);
    for i = 1:5
        x_rho(:,i) = x_boundary(s,cmpd,(2*i-1)/10)./s.EL;
    end
    mean_x_rho = mean(x_rho);
    std_x_rho = std(x_rho);
    % Add error bars in position at normalized intensities of 0.1, 0.3, 0.5,
    % 0.7, 0.9
    errorbar(mean_x_rho,(2*(1:5)-1)/10,std_x_rho,'horizontal','k.','LineWidth',2);
    xlim([0 1])
    ylim([0 1])
end

function stat_plot(s_arr,x_vec)
% Plots the mean position of the Hb boundary for a number of simulation
% trials
% Inputs
    % s_arr - Array of structs containing data from different simulation
    % trials
    % x_vec - Vector of positions to be plotted on the x-axis
    
    x_Hb = zeros(s_arr(1).N,length(s_arr));
    for i = 1:length(s_arr)
        x_Hb(:,i) = x_boundary(s_arr(i),5,0.5)./s_arr(i).EL;
    end
    stats_x_Hb = [mean(x_Hb); std(x_Hb)];
    errorbar(x_vec,stats_x_Hb(1,:),stats_x_Hb(2,:))
end
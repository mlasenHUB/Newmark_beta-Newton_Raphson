%% Plot figures for the Lumped model
close all
clear all
xt_min = 0.09;
xt_max = 0.14;
yd_max = 0.012;
dofs = 3;%plate nodes
dp = 4; %punch
dmm = 5; %moving mass
%% fig 100
%load('MDF_with_plate_k_const_0  1  0.mat')
%load('E:\From Time Integration\MDF_with_plate_k_const_and_diff_k_plate0  1  0.mat')
%load('E:\From Time Integration\MDF_with_plate_adjusted_kt_to simulations_C6_and_diff_k_plate0  1  0.mat')
%load('E:\From Time Integration\MDF_with_plate_extracted_kt_from_simulations_C6_and_diff_k_plate0  1  0.mat', 't', 'u_x', 'f_f')
load('E:\From Time Integration\test_0  1  0.mat', 't', 'u_x', 'f_f')

figure(100)

plot(t, 1e3*u_x(4:5, :), 'LineWidth', 5)
hold on
plot(t, 1e3*1e2*u_x(1:3, :), 'LineWidth', 5)
legend('DOF 4', 'DOF 5','DOF 1 x100', 'DOF 2 x100','DOF 3 x100')
 
grid on
xlabel('Time [s]')
ylabel('Displacement[mm]')
ylim([-yd_max yd_max])
xlim([xt_min, xt_max])
set(gca,'FontSize',40)
set(gcf, 'Position', get(0, 'Screensize'));

figure(103)
in_ind = length(u_x(1, :))-4.3*1e5;
u_mean_plate = mean(u_x(1:dofs, in_ind:end),1);
rel_dis = u_x(dp, in_ind:end)-u_mean_plate;
plot(rel_dis*1e3, f_f(dmm-1, in_ind:end), 'DisplayName', '010','LineWidth', 5)
hold on

figure(104)
plot(t, f_f(dmm-1, :),'DisplayName', '010', 'LineWidth', 5)
hold on

%% fig 101
%load('MDF_with_plate_k_const_1  0  1.mat')
%load('E:\From Time Integration\MDF_with_plate_k_const_and_diff_k_plate1  0  1.mat')
%load('E:\From Time Integration\MDF_with_plate_adjusted_kt_to simulations_C6_and_diff_k_plate1  0  1.mat')
%load('E:\From Time Integration\MDF_with_plate_extracted_kt_from_simulations_C6_and_diff_k_plate1  0  1.mat','t', 'u_x', 'f_f')
load('E:\From Time Integration\test_1  0  1.mat', 't', 'u_x', 'f_f')

figure(101)
plot(t, 1e3*u_x(4:5, :), 'LineWidth', 5)
hold on
plot(t, 1e3*1e2*u_x(1:3, :), 'LineWidth', 5)
legend('DOF 4', 'DOF 5','DOF 1 x100', 'DOF 2 x100','DOF 3 x100')
 
grid on
xlabel('Time [s]')
ylabel('Displacement[mm]')
ylim([-yd_max yd_max])
xlim([xt_min, xt_max])
set(gca,'FontSize',40)
set(gcf, 'Position', get(0, 'Screensize'));
figure(103)
in_ind = length(u_x(1, :))-4.2*1e5;
u_mean_plate = mean(u_x(1:dofs, in_ind:end),1);
rel_dis = u_x(dp, in_ind:end)-u_mean_plate;
plot(rel_dis*1e3, f_f(dmm-1, in_ind:end), 'DisplayName', '101','LineWidth', 5)

figure(104)
plot(t, f_f(dmm-1, :),'DisplayName', '101', 'LineWidth', 5)

%% fig 102
%load('MDF_with_plate_k_const_1  1  1.mat')
%load('E:\From Time Integration\MDF_with_plate_k_const_and_diff_k_plate1  1  1.mat')
%load('E:\From Time Integration\MDF_with_plate_adjusted_kt_to simulations_C6_and_diff_k_plate1  1  1.mat')
%load('E:\From Time Integration\MDF_with_plate_extracted_kt_from_simulations_C6_and_diff_k_plate1  1  1.mat', 't', 'u_x', 'f_f')
load('E:\From Time Integration\test_1  1  1.mat', 't', 'u_x', 'f_f')

figure(102)
plot(t, 1e3*u_x(4:5, :), 'LineWidth', 5)
hold on
plot(t, 1e3*1e2*u_x(1:3, :), 'LineWidth', 5)
legend('DOF 4', 'DOF 5','DOF 1 x100', 'DOF 2 x100','DOF 3 x100')
 
grid on
xlabel('Time [s]')
ylabel('Displacement[mm]')
ylim([-yd_max yd_max])
xlim([xt_min, xt_max])
set(gca,'FontSize',40)
set(gcf, 'Position', get(0, 'Screensize'));


figure(103)
in_ind = length(u_x(1, :))-5*1e5;
u_mean_plate = mean(u_x(1:dofs, in_ind:end),1);
rel_dis = u_x(dp, in_ind:end)-u_mean_plate;
plot(rel_dis*1e3, f_f(dmm-1, in_ind:end), 'DisplayName', '111','LineWidth', 5)

figure(104)
plot(t, f_f(dmm-1, :),'DisplayName', '111', 'LineWidth', 5)

%% fig 103
%All Hysteresis Loops 010/101/111
xli = 15e-3;%[mm]
figure(103)
hold on
load('E:\From Time Integration\MDF_without_plate_60N_tang_Load.mat')
in_ind= length(u_x(1, :)) - 1*1e5;
plot(u_x(1, in_ind:end)*1e3, f_f(1, in_ind:end), 'DisplayName', 'No plate','LineWidth', 5)
legend show

xlabel('Relative Displacement [mm]')
ylabel('Tangential Force [N]')
xlim([-xli xli])
set(gca,'FontSize',40)
grid on
hold off
set(gcf, 'Position', get(0, 'Screensize'));

%% fig 4
figure(104)
hold on

load('E:\From Time Integration\MDF_without_plate_60N_tang_Load.mat', 't', 'p_p', 'f_f')
plot(t, f_f(dmm-1, :),'DisplayName', 'No plate', 'LineWidth', 5)
ylabel('Tangential Force [N]')
yyaxis right
plot(t, p_p, 'DisplayName', 'Excitation', 'LineWidth', 5)
ylabel('Excitation Force [N]')
legend show
xlabel('Time [s]')

xlim([0 xt_max])
ylim([-70 70])
set(gca,'FontSize',40)
grid on
hold off
set(gcf, 'Position', get(0, 'Screensize'));

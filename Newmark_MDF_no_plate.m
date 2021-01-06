%%% Newmark-beta for MDOF with friction nonlinearity
%%% Using Newton-Raphson method with constant Jacobian
%%% Different patterns incorporated  (010,101,111, in gral)
%%% Patterns modifiable in time (2 of them)
%%% With Convergence criteria
%%% 3 dofs Plate - 2dofs moving mass+plate+punch

clear all
close all

%% System properties

dofs = 0; %no plate
dp = dofs + 1; %punch degree of freedom
dmm = dofs + 2; %add the punch

m_mm = 21.2;%moving mass
m_p = 0.01163+1.2; %mass of the punch+ moving arm


M = zeros(dmm);
M(dmm-1, dmm-1) = m_p;
M(dmm, dmm) = m_mm;

k_p = 280e6;%moving arm to moving mass 273
k_mm = 1.35e6;%moving mass to ground 1.35

K = make_K(k_p, k_mm, 2);

w_n_p = sqrt(k_p/m_p); %nat freq of the punch
w_n_mm = sqrt(k_mm/m_mm); %nat freq of the punch
% c_p = chi_p*2*m_p*w_n_p; % moving arm to moving mass
% c_mm = chi_mm*2*m_mm*w_n_mm; % moving mass to ground

c_p = 0.0151; chi_p = c_p/(2*m_p*w_n_p);%0.0151
c_mm = 0.006; chi_mm = c_mm/(2*m_mm*w_n_mm);%c_mm = 0.006
C = make_C(c_p, c_mm, 2);

%% Contact Stiffness
k_t_exp = 50e6;%50e6experimental k_t steel/steel
k_t = 1*k_t_exp; %need to change if the first pattern is changed
K_t = [ k_t 0; 0 0];

%% Eigen - Stuck Punch
K_st = K + K_t;
[V_st,D_st] = eig(K_st,M);%V:eigen vectors, D:eigen values.
d_st=diag(D_st); 
W_n_st = sqrt(d_st);%[2pi/s]
f_n_st = W_n_st/2/pi; %[Hz]

e_st = dmm;
r_st = dp;
f_st = zeros(1, dmm);
f_st(e_st) = 1;
% [ Hhast, Hhtast, phst, Vst, Dst, fdomst ] = LinearFrequencyResponse( M,K_st, C, f_st, r_st, k_t);

%% Time parameters
T_n_s = 1/W_n_st(dmm); % [Hz], s: entire System, stuck

dt = T_n_s/10;
t_end = 20000*T_n_s; %Number of periods
t = 0:dt:t_end; % total length of time

n = length(t)-1; % number of time steps
t = t(1:n);
%% Tangential Load
P0 = 60; %micro N
w_ext = 100;
p_p = P0*sin(2*pi*w_ext*t);
P = zeros(dmm, n);
P(dmm, :) = p_p; %Excitatioin in the moving mass
%% Newmark -beta cases
method=1; %Const.ave.accel., else linear
if dt/T_n_s<=0.551 && method ~= 1 % Use linear accel. method - closest to theoretical
    gamma=0.5;
    beta=1/6;
end
if dt/T_n_s<=0.318 && method==1 % Use constant avg. accel. method - unconditionally stable (Example 5.5)
    gamma=0.5;
    beta=0.25;
end
%% Initial conditions
u_x = zeros(dmm, n);
v_x = zeros(dmm, n);
a_x = zeros(dmm, n);
%% Newmark constants
A1_x = M*(1/(beta*dt*dt)) + C*(gamma/(beta*dt));
A2_x = M*(1/(beta*dt))+ C*(gamma/beta-1);
A3_x = M*((0.5/beta)-1) + C*dt*(gamma/(2*beta)-1);

p_hat_x = zeros(dmm,n);
R = zeros(dmm, n);
R_tol = 10^-3;
%% Normal Load - Friction

mu = 0.67;
w_N0 = 0*w_ext;% 0*w_ext means constant N0
max_N0 = 40;
min_N0 = 15;
Amp_pp = max_N0-min_N0; %Amplitude peak to peak

N0 = Amp_pp/2*(cos(2*pi*w_N0*t)+1)+min_N0; %total normal force
f_f = zeros(dmm, n); %restoring force
%% Run cases
number_case = 1;
switch  number_case
    
    case 1 %
        close all
        %% Initial Conditions
        a_x(:,1) = M\( P(:,1) - C*v_x(:,1) - K*u_x(:,1) - f_f(:,1));
%         conv =0; %converged =1
        K_t_old = K_t;%there is no previous iter for the first increment
        
        
        for i=1:n-1 %time
            
%             if conv==1
%                 break%stop iteration if converged to steady
%             end
            p_hat_x(:,i+1) = P(:,i+1)+ A1_x*u_x(:,i)+ A2_x*v_x(:,i)+...
                A3_x*a_x(:,i);
            
            %% Residual first iteration
            u_x(:,i+1) = u_x(:,i);
            f_f(:,i+1) = f_f(:,i);
            
            R(:,i+1) = p_hat_x(:,i+1)- K*u_x(:,i+1)- f_f(:,i+1)-...
                A1_x*u_x(:,i+1);
            
            %% Update pattern, stiffness matrix,normal load and friction limits
            
            fr_limit = mu*N0(i);
            %1st increment is stuck
            K_t = [ k_t 0; 0 0];
            
            if K_t ~=K_t_old
                K_t =K_t_old; %if any dof was previously sliding, keep it like that for the next increment
            end
            
            %% Newton -Raphson
                       
            while norm(R(:,i+1))> R_tol
           
                J = A1_x + K + K_t;
                delta_x = J\R(:,i+1);
                u_x(:,i+1) = u_x(:,i+1) + delta_x;
                
                f_f(:, i+1) = f_f(:, i+1) + K_t*delta_x;
%
                %SLIP<->STICK Transitions for each DOF
                %stick->slip
                if k_t ~=0 && abs(f_f(1,i+1)) > fr_limit
                    %positive
                    if f_f(1,i+1) > fr_limit
                        f_f(1,i+1)= fr_limit;
                    end
                    %negative
                    if f_f(1,i+1) < -fr_limit
                        f_f(1,i+1)= -fr_limit;
                    end
                    %update K_t
                    K_t(1, 1) = 0;
                end
                %slip->stick
                if k_t ~=0 && abs(f_f(1,i+1)) <= fr_limit
                    K_t(1, 1) = k_t;
                end
                K_t_old = K_t;%for the new increment this LAST iteration will be the old case
                 
                %update residual
                R(:,i+1) = p_hat_x(:,i+1)- K*u_x(:,i+1)- f_f(:,i+1)- A1_x*u_x(:,i+1);
            end
            % update final values
        v_x(:,i+1) = gamma*(u_x(:,i+1)-u_x(:,i))/(beta*dt) + v_x(:,i)*(1-gamma/beta) + dt*a_x(:,i)*(1-gamma/(2*beta));
        a_x(:,i+1) = (u_x(:,i+1)-u_x(:,i))/(beta*dt*dt) - v_x(:,i)/(beta*dt) -((0.5/beta) -1)*a_x(:,i);
        
        end
        
        
end

%% FIGURE 4
        delta = 2500;
        init = length(u_x(1, :)) - delta;
        %init =1;
        figure(4)
        subplot(1, 4, 1)
        plot(t(init:end), u_x(:,init:end))
        legend('u_x_1', 'u_x_2')
        
        subplot(1, 4, 2)
        plot(t(init:end), u_x(dmm-1, init:end))
        hold on
        xlabel('Time [s]')
        hold off
        yyaxis right
        plot(t(init:end), f_f(dmm-1, init:end))
        grid on
        ylabel('Force [N]')
        legend('u_x_1', 'Tangential Force in the Punch')
        
        subplot(1, 4, 3)
        plot(t(init:end), u_x(dmm-1,init:end))
        yyaxis right
        plot(t(init:end), v_x(dmm-1,init:end))
        legend('u_x_1', 'v_x_1')
        
        subplot(1, 4, 4)
        plot(t(init:end), N0(init:end))
        hold on 
        plot(t(init:end), p_p(init:end))
        plot(t(init:end), f_f(2, init:end))
        grid on
        xlabel('Time [s]')
        ylabel('Force [N]')
        legend('Total Normal Load [N]', 'Tangential Load [N]', 'Friction force on the punch [N]')
         %% FIGURE 6
         
         figure(6)
         plot(u_x(1, init:end), f_f(1, init:end), 'o','DisplayName', 'Punch', 'LineWidth', 5)
         xlabel('Displacement [m]')
         ylabel('Tangential Force [N]') 
         set(gca,'FontSize',25)
         %% Figure 7
         figure(7)
         plot(t(init:end), u_x(:,init:end),'LineWidth', 5)
         hold on
         yyaxis right
         plot(t(init:end), f_f(1, init:end), 'LineWidth', 5)
         hold off
         legend('u_x_1', 'u_x_2', 'f_f_1')
         set(gca,'FontSize',25)
save('E:\From Time Integration\MDF_without_plate_60N_tang_Load')


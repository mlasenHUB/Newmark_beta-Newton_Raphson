%%% Newmark-beta for MDOF with friction nonlinearity
%%% Using Newton-Raphson method with constant Jacobian
%%% Different patterns incorporated  (010,101,111, in gral)
%%% Patterns modifiable in time (2 of them)
%%% With Convergence criteria
%%% 3 dofs Plate - 2dofs moving mass+plate+punch

clear all
close all
%% System properties

dofs = 3;
dp = dofs + 1; %punch degree of freedom
dmm = dofs + 2; %add the punch
%k_cons = 1; %true k constant for al elements, 0: false
%if false
%k_1 = 1.5; k_2 =1.2; k_3 =1; %multipliers: k_1 if ONE element in contact

m_plate = 0.00032;
m_lat = m_plate*(8.5/19);
m_center = m_plate*(2/19);
m_mm = 21.2;%moving mass
m_p = 0.01163+1.2; %mass of the punch+ moving arm

m_lat = m_plate/3; m_center=m_lat;
M = diag([m_lat m_center m_lat]);
M_1D = diag([m_p m_mm]);
M(dp:dmm,dp:dmm) = M_1D;

kL = 63.2e6;%steel plate stiffnes kL= EA/L, L=19mm
k1 = kL*(19/7.5);
k2 = kL*(19/2);
%k1 =4*kL; k2=k1;
k_p = 280e6;%moving arm to moving mass 273
k_mm = 1.35e6;%moving mass to ground 1.35e6


K = make_K(k1, k2, 3);
K_1D = make_K(k_p, k_mm, 2);
K(dp:dmm, dp:dmm) = K_1D;


k = min(k1, k2);
m = max(m_lat, m_center);
w_n = sqrt(k/m); %nat freq of the plate
w_n_p = sqrt(k_p/m_p); %nat freq of the punch
w_n_mm = sqrt(k_mm/m_mm); %nat freq of the punch
% chi = 2; %Damping ratio plate
% chi_p = 1; %0.6 from alfredo FRF
% chi_mm = 1;
% c = chi*2*m*w_n;
% c_p = chi_p*2*m_p*w_n_p; % moving arm to moving mass
% c_mm = chi_mm*2*m_mm*w_n_mm; % moving mass to ground

c_p = 0.0151; chi_p = c_p/(2*m_p*w_n_p);
c_mm = 0.006; chi_mm = c_mm/(2*m_mm*w_n_mm);
c =100; chi = c/(2*m*w_n);

c1=c ; c2=c;
C = make_C(c1, c2, 3);
C_1D = make_C(c_p, c_mm, 2);
C(dp:dmm, dp:dmm) = C_1D;
%% Pattern initialisation
pattern1 = [1 1 1];
pattern2 = pattern1; %don't change

pattern = pattern1;
%% Contact Stiffness 
%proportional to contact area total area

% tot_area_010 = 1.168; %from simulation t_82
% tot_area_101 = 4.545; %from simulation t_83_new
% tot_area_111 = 8.295; %from simulation t_84_new

%Extract average areas in time from simulaton with clear pattern(need to
%see voltage)
%NOTE THAT STILL DOES NOT EXPLICITELY TAKE INTO ACCOUNT THE EFFECT OF THE
%DIFFERENT NORMAL LOAD PER CONTACT POINT.
tot_area_010 = 3; %from simulation C6_010_80
tot_area_101 = 6.7; %from simulation C6_101_50
tot_area_111 = 9.6; %from simulation C6_101_50
k_t_exp = 50e6;%50e6experimental N/m k_t steel/steel for 1mm2 contact area 

% k_ts per element
if pattern == [0 1 0]
    %k_t = k_t_exp*tot_area_010/sum(pattern);
    k_t = 7.313*1e6;%from C6_010_80 Hysteresis Loop
elseif pattern == [1 0 1]
    %k_t = k_t_exp*tot_area_101/sum(pattern);
    k_t = 4.3198*1e6;%from C6_101_80 Hysteresis Loop
elseif pattern == [1 1 1]
    %k_t = k_t_exp*tot_area_111/sum(pattern);
    k_t = 3.3234*1e6;%from C6_111_80 Hysteresis Loop
end

k_t_v = k_t*pattern';
k_t_v_h = k_t*[pattern 0];
k_t_m = k_t*diag(pattern);

K_t= [ k_t_m -k_t_v];
K_t = [K_t; -k_t_v_h ];
K_t(dmm-1, dmm-1) = sum(pattern)*k_t;
K_t(dmm, dmm) = 0;

%% Eigen - Stuck Punch 
K_st = K + K_t;
[V_st,D_st] = eig(K_st,M);%V:eigen vectors, D:eigen values.
d_st=diag(D_st); chi_m = c/2/d_st(1); %damping coeficient times mass
W_n_st = sqrt(d_st);

e_st = 5;
r_st = 4;
f_st = zeros(1, dmm);
f_st(e_st) = 1;
% [ Hhast, Hhtast, phst, Vst, Dst, fdomst ] = LinearFrequencyResponse( M,K_st, C, f_st, r_st, k_t);
%% Time parameters
T_n_s = 1/W_n_st(dmm); % [Hz], s: entire System, coupled
T_n = 1/w_n;%plate only
T_n_p = 1/w_n_p;%punch only
dt = T_n_s/10;

t_end = 600000*T_n_s; %Number of periods
t = 0:dt:t_end; % total length of time

n = length(t)-1; % number of time steps
t = t(1:n);
%% Tangential Load
P0 = 60;
w_ext = 100;
% unitstep = P*(t>=0);
% unitstep2 = P*(t>=1*dt);
p_p = P0*sin(2*pi*w_ext*t);
P = zeros(dmm, n);
P(dmm, :) = p_p; %Excitatioin in the moving mass
%P = repmat(p_x, dofs, 1); %repeat same excitation for all dofs
%% Pattern - Activation 

f =4; %Multiple of the excitation freq.
w_ext_pattern = f*w_ext;%frequency of alternation between pattern
%alternate patterns
sqwave = 0.5*square(2*pi*w_ext_pattern*t);
sqwave =sqwave +0.5;
sqwave2 = 0.5*square(2*pi*w_ext_pattern*(t)- pi);
sqwave2 = sqwave2+0.5;
%generate a single vector with the alternating patterns
pattern_v = pattern1'*sqwave + pattern2'*sqwave2;
pattern_v = pattern_v'; %patterns in time are the rows

cont_s = sum(pattern_v,2);% number of elements in contact at any time
% k_t_V = zeros(1, length(cont_s));% corresponding values of the tang stiffnesses
% 
% if k_cons==0
%     for i = 1:length(cont_s)
%         if cont_s(i) ==1
%             k_t_V(i)= k_1;%1.5
%         elseif cont_s(i) ==2
%             k_t_V(i)= k_2;%1.2
%         else
%             k_t_V(i)= k_3;
%         end
%     end
% end
%if k_cons==1
k_t_V = ones(1, length(cont_s)); %vector containing kt at each time for every element in contact, for any pattern
%end

%% Contact Stiffness - Activation
%find the value of the tangential stiffness of the second pattern(if exist)
aux_ind = find(k_t_V~=k_t_V(1),1);
k_t2 = k_t_V(aux_ind);

if isempty(k_t2)
    k_t2=k_t;% set it as the same as the first if there is none
end

k_t_v2 = k_t2*pattern2';
k_t_v_h2 = k_t2*[pattern2 0];
k_t_m2 = k_t2*diag(pattern2);

K_t2= [ k_t_m2 -k_t_v2];
K_t2 = [K_t2; -k_t_v_h2 ];
K_t2(dmm-1, dmm-1) = sum(pattern2)*k_t2;
K_t2(dmm, dmm) = 0; %K_t for the second pattern

K_t1 = K_t; % K_t for the first pattern
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
R = zeros(dmm, n); %residual initialization
R_tol = 10^-3;
%% Normal Load - Friction

mu = 0.67;
w_N0 = 0*w_ext;% 0*w_ext means constant N0
max_N0 = 40;
min_N0 = 15;%for variable Normal Load
Amp_pp = max_N0-min_N0; %Amplitude peak to peak

N0_T = Amp_pp/2*(cos(2*pi*w_N0*t)+1)+min_N0; %total normal force
f_f = zeros(dmm, n); %friciton force
mu_v = mu*ones(1, dofs);

N0= N0_T./cont_s'; %normal load per element
k_t_V = k_t_V*k_t; %k_ts defined as proportoinal to area or not
%% Run cases
number_case = 1;
switch  number_case
    
    case 1 %
        close all
        %% Initial Conditions
        a_x(:,1) = M\( P(:,1) - C*(v_x(:,1)) - K*(u_x(:,1)) - f_f(:,1));
        K_t_old = K_t1;%set the first tang stiffness mat
        for i=1:n-1 %time

            p_hat_x(:,i+1) = P(:,i+1)+ A1_x*u_x(:,i)+ A2_x*v_x(:,i)+...
                A3_x*a_x(:,i);
            
            %% Residual first iteration
            u_x(:,i+1) = u_x(:,i);%Displacement
            f_f(:,i+1) = f_f(:,i);%Friction
            R(:,i+1) = p_hat_x(:,i+1)- K*u_x(:,i+1)- f_f(:,i+1)-...
                A1_x*u_x(:,i+1); 
            %% Update pattern, stiffness matrix,normal load and friction limits
            pattern = pattern_v(i, :);%what is the current pattern
            k_tt = k_t_V(i);% local contact stiffness of any of the elements in the pattern, for any pattern
            fr_limit_v = mu_v*N0(i).*pattern;
            
            %1st increment is stuck
%             if cont_s(i) == sum(pattern1)
%                 K_t = K_t1;
%             end
%             if cont_s(i) == sum(pattern2)
%                 K_t = K_t2;
%             end
%             
            if K_t ~=K_t_old
                K_t =K_t_old; %if any dof was previously sliding, keep it like that for the next increment
            end
         %% Newton -Raphson
            while norm(R(:,i+1)) > R_tol
                
                J = A1_x + K + K_t;
                delta_x = J\R(:,i+1);
                u_x(:,i+1) = u_x(:,i+1) + delta_x;
                %Restoring force in the plate Dofs
%                 f_s(1:dofs,i+1) = f_s(1:dofs,i+1) + K_t(1:dofs, 1:dofs)*...
%                     (u_x(1:dofs,i+1)-u_x(1:dofs,i)-u_x(dmm-1, i+1)); % within the iteration
%                 f_s(1:dofs,i+1) = f_s(1:dofs,i+1) + K_t(1:dofs, 1:dofs)*...
%                     (u_x(1:dofs,i+1)-u_x(1:dofs,i)-u_x(dmm-1, i+1)); % within the iteration
%                 
                f_f(:, i+1) = f_f(:, i+1) + K_t*delta_x;
%
                %SLIP<->STICK Transitions for each DOF
                for ii = 1:dofs
                    fr_limit = fr_limit_v(ii);%for each dofs
                    
                    %stick->slip
                    if k_tt ~=0 && abs(f_f(ii,i+1)) > fr_limit
                        %positive
                        if f_f(ii,i+1) > fr_limit
                            f_f(ii,i+1)= fr_limit;                            
                        end
                        %negative
                        if f_f(ii,i+1) < -fr_limit  
                            f_f(ii,i+1)= -fr_limit;                            
                        end 
                        %update K_t for slipping condition                       
                        K_t(ii, ii) = 0;
                        K_t(ii, dp) = 0; 
                        K_t(dp, ii) = 0;                        
                    end
                    %slip->stick
                    if k_tt ~=0 && abs(f_f(ii,i+1)) <= fr_limit 
                        if pattern(ii) ==1%pattern =1 =>jenkins element actually exist
                        K_t(ii, ii) = k_tt;
                        K_t(ii, dp) = -k_tt;
                        K_t(dp, ii) = -k_tt;
                        end
                    end
                end
                K_t(dp, dp) = -sum(K_t(dp, :));
                K_t_old = K_t;%for the new increment this NEW iteration will be the old case
                % Friction in Punch is the total sum of all others
                f_f(dmm-1, i+1) = -sum(f_f(1:dofs, i+1));
                %update residual
                R(:,i+1) = p_hat_x(:,i+1)- K*u_x(:,i+1)- f_f(:,i+1)- A1_x*u_x(:,i+1);
            end
            
            % update final values after convergence criteria (Residual) is
            % met
            v_x(:,i+1) = gamma*(u_x(:,i+1)-u_x(:,i))/(beta*dt) + v_x(:,i)*(1-gamma/beta) + dt*a_x(:,i)*(1-gamma/(2*beta));
            a_x(:,i+1) = (u_x(:,i+1)-u_x(:,i))/(beta*dt*dt) - v_x(:,i)/(beta*dt) -((0.5/beta) -1)*a_x(:,i);
            %% Convergence Check
            
%             n_last = 6;
%             c_c = 1; % percent for the convergence criteria
%             inc_per = 40000;%increment periodicity
            %if i==inc_ind*inc_per
%                 [ave_curr, conv, dif] = convergence(u_x(4, 1:i), t(1:i), n_last, c_c);
% %                 conv, dif, i
%                 inc_ind = inc_ind + 1;
%                 if conv==1
%                     break
%                 end
%             end
            
        end 
        set(gca,'FontSize',25)
        %% Figure 1
        figure (1)
        suptitle (['Patterns: ', num2str(pattern1), ' >> ',num2str(pattern2)])
        subplot(2, 3, 1)
        for i =1:dofs
            txt = ['v_x_', num2str(i)];
            plot(t, v_x(i,:), 'DisplayName', txt)
            hold on
        end
        hold off
        legend show
        xlabel('Time [s]')
        ylabel('Velocity  [m/s]')
        
        subplot(2, 3, 2)
        for i =1:dofs
            txt = ['u_x_', num2str(i)];
            plot(t, u_x(i,:), 'DisplayName', txt)
            hold on
        end
        hold off
        legend show
        xlabel('Time [s]')
        ylabel('Displacement [m]')
        
        subplot(2, 3, 3)
        for i =1:dofs
            txt = ['f_s_', num2str(i)];
            plot(t, f_f(i,:), 'DisplayName', txt)
            hold on
        end
        
        txt = [txt, 'Normal Load'];
        legend show
        xlabel('Time [s]')
        ylabel('Tangential Force [N]')
        yyaxis right
        plot(t, N0_T)
        ylabel('Normal Load [N]')
        hold off
        
        subplot(2, 3, 4)        
        plot(t, v_x(dmm,:))
        legend('Punch')
        xlabel('Time [s]')
        ylabel('Velocity  [m/s]')
        
        subplot(2, 3, 5)        
        plot(t, u_x(dmm,:))
        legend('Punch')
        xlabel('Time [s]')
        ylabel('Displacement [m]')
        
        subplot(2, 3, 6)        
        plot(t, f_f(dmm-1,:))
        ylabel('Tangential Force [N]')
        
        yyaxis right
        hold on
        plot(t, P(dmm-1,:))
        plot(t, sqwave)
        plot(t, sqwave2)
        legend('Punch Contact force [N]', 'Punch Excitation force [N]', 'Pattern 1', 'Pattern 2')
        xlabel('Time [s]')
        hold off
        
        %% Figure 2
        figure(2)
        plot(t, u_x(1:3,:), 'LineWidth', 5) 
        xlim([0, 0.05])
        legend(' DOF 1', 'DOF 2', 'DOF 3')
        
        %% FIGURE 4
        figure(4)
        subplot(1, 4, 1)
        plot(t, u_x)
        legend('u_x_1', 'u_x_2','u_x_3','u_x_p', 'u_mov_m')
        
        subplot(1, 4, 2)
        %plot(t, u_x(1:3, :))
        plot(t, u_x(dmm-1, :))
        %hold on
        %plot(t, v_x(dofsp-1, :))
        xlabel('Time [s]')
        %ylabel('Displacement [m]')
        %hold off
        yyaxis right
        %plot(t, u_x(4, :))        
        plot(t, f_f(dmm-1, :))
        grid on
        ylabel('Force [N]')
        legend('u_x_p', 'Tangential Force in the Punch')
        
        subplot(1, 4, 3)
        plot(t, u_x(dmm-1,:))
        yyaxis right
        plot(t, v_x(dmm-1,:))
        legend('u_x_p', 'v_x_p')
        
        subplot(1, 4, 4)
        plot(t, N0_T)
        hold on 
        plot(t, p_p)
        plot(t, f_f(4, :))
        grid on
        xlabel('Time [s]')
        ylabel('Force [N]')
        legend('Total Normal Load [N]', 'Tangential Load [N]', 'Friction force on the punch [N]')
        hold off
        suptitle (['Patterns: ', num2str(pattern1), ' >> ',num2str(pattern2), '; Normal Load/Maximum Normal Load = ', num2str(max_N0), '[N]; Excitation Amplitude = ', num2str(P0), '[N]'])
             
        
%         figure(5)
%         for i =1:dofsp
%             subplot(1, dofsp, i)
%             txt = ['Hysteresis-DOF_', num2str(i)];
%             plot(u_x(i, :), -1*f_s(i, :), 'DisplayName', txt)
%             legend show
%             xlabel('Displacement [m]')
%             ylabel('Tangential Force [N]')
%         end
        %% FIGURE 6
        figure(6)
        rel_dis_plate = u_x(1:3, :)-u_x(4,:);
        in_ind = 1;%n-round(n/10);
        xli = 15e-3;%[mm]
        for i =1:dofs
            subplot(dmm-1,1, i)
            txt1 = ['Hysteresis-DOF-RELATIVE', num2str(i)];
            txt2 = ['Hysteresis-DOF-SELF', num2str(i)];
            plot((u_x(i, in_ind:end)-u_x(dp,in_ind:end))*1e3, -1*f_f(i, in_ind:end), 'DisplayName', txt1, 'LineWidth', 3)
            hold on
            plot(u_x(i, in_ind:end)*1e3, -1*f_f(i, in_ind:end), 'DisplayName', txt2, 'LineWidth', 3)
            legend show
            xlabel('Displacement [mm]')
            ylabel('Tangential Force [N]')
            xlim([-xli xli])
            set(gca,'FontSize',25)
            
        end
        hold off
        subplot(dmm-1,1,  dmm-1)
        %plot(u_x(dp, in_ind:end)-u_x(2,in_ind:end), f_s(dmm-1, in_ind:end))
        u_mean_plate = mean(u_x(1:dofs, in_ind:end),1);
        rel_dis = u_x(dp, in_ind:end)-u_mean_plate;
        plot(rel_dis*1e3, -f_f(dmm-1, in_ind:end), 'LineWidth', 3)
        hold on
        plot(u_x(dp, in_ind:end)*1e3, f_f(dmm-1, in_ind:end), 'LineWidth', 3)
        hold off
        legend('Hys_{Average}','Hys_{Absolute}' )
        xlabel('Relative Displacement [mm]')
        ylabel('Tangential Force [N]')
        xlim([-xli xli])
        set(gca,'FontSize',25)
        
        %% figure 7 
        figure(7)
        
        yyaxis left
        plot(t, rel_dis*1e3, 'DisplayName', 'Relative Displacement', 'LineWidth', 5)
        ylabel('Displacement [mm]')
        ylim([-0.012 0.012])
        ax = gca;
        yyaxis right
        plot(t, p_p,'DisplayName', 'Excitation Force','LineWidth', 5)
        hold on
        plot(t, f_f(4, :), 'DisplayName', 'Friction Force','LineWidth', 5)
        hold off
        ax = gca;
        ax.YAxis(2).Color = [0.9290 0.6940 0.1250];
        ylabel('Force [N]')
        ylim([-1.1*P0 1.1*P0])
        xlim([0 0.1])
        xlabel('Time [s]')
        legend show
        grid on
        set(gca,'FontSize',25)
        
        %% FIGURE 8
        figure(8)
        
        plot(t, u_x(4:5, :))
        hold on
        plot(t, 1e2*u_x(1:3, :))
        legend('u_x_p', 'u_mov_m','u_x_1', 'u_x_2','u_x_3')
        xlabel('Time [s]')
        legend show
        grid on
        set(gca,'FontSize',25)
        hold off
%         %% FIGURE 8
%         figure(8)
%         for i =1:dmm
%             txt = ['Residual-DOF_', num2str(i)];
%             plot(R(i, :), 'DisplayName', txt)
%             hold on
%             legend show
%         end
        
        
        % Frequency results
%         figure(9)
%         for i=1:dofs
%             
%             txt = ['Mode_' , num2str(i)];
%             semilogy(fdom,abs(Hha(:,i)), 'DisplayName', txt)
%             hold on
%         end
%         grid on
%         
%         
%         semilogy(fdom,abs(Hhta), 'DisplayName', ['FRF_e_', num2str(e), '_r_', num2str(r)])
%         legend show
%         title('3 DOFs Response - No Punch')
%         xlabel('Frequency [Hz]')
%         
%         for i=1:dofs
%             
%             txt = ['Stuck Mode_' , num2str(i)];
%             semilogy(fdoms,abs(Hhas(:,i)), 'DisplayName', txt)
%             
%         end
%         grid on
%         
%         semilogy(fdoms,abs(Hhtas), 'DisplayName', ['Stuck FRF_e_', num2str(e), '_r_', num2str(r)])
%         legend show
%         
%         hold off
        
%         figure(10)
%         for i=1:dofsp
%             
%             txt = ['Mode_' , num2str(i)];
%             semilogy(fdomst,abs(Hhast(:,i)), 'DisplayName', txt)
%             hold on
%         end
%         grid on
%         
%         
%         semilogy(fdomst,abs(Hhtast), 'DisplayName', ['FRF_e_', num2str(e_st), '_r_', num2str(r_st)])
%         legend show
%         title('4 DOFs Response - With Punch')
%         xlabel('Frequency [Hz]')
        
%         figure(11)
%         findpeaks(u_x(4, :))
          
end
% % file = ['E:\From Time Integration\MDF_with_plate_extracted_kt_from_simulations_C6_and_diff_k_plate',...
%     num2str(pattern)];
file = ['E:\From Time Integration\test_',...
    num2str(pattern)];
save(file)

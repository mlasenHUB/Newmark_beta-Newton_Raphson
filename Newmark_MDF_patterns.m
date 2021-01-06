%%% Newmark-beta for MDOF with friction nonlinearity
%%%Using Newton-Raphson method with constant Jacobian
%%% Different patterns incorporated

clear all
close all


%% Establish system properties

dofs = 3;
dofsp = dofs + 1; %add the punch

m = 1;
m_p = 3; %mas of the punch
M_d = m*ones(1, dofsp);
M = diag(M_d);
M(dofsp, dofsp) = m_p;


k = 100;
k_p = 100;
K = diag(2*k*ones(1, dofsp)) - diag(k*ones(1,dofsp-1),1) - diag(k*ones(1, dofsp-1),-1);
K(dofsp, dofsp) = k_p;
K(dofsp, dofs) = 0;
K(dofs, dofsp) = 0;


c = 0.1;
c_p = 0.3;
C = diag(2*c*ones(1, dofsp)) - diag(c*ones(1,dofsp-1),1) - diag(c*ones(1, dofsp-1),-1);
C(dofsp, dofsp) = c_p;
C(dofsp, dofs) = 0;
C(dofs, dofsp) = 0;

% %% Linear Response of Membrane Only ( NO PUNCH)
% e = 1; %Excitation point 
% r = 2; %Response point
% f = zeros(1, dofs);
% % Excitation force (only for linear analysis)
% f(e) = 1;
% f(3) = -1;
% 
% %Free Response
% [ Hha, Hhta, ph, V, D, fdom ] = LinearFrequencyResponse( M(1:dofs, 1:dofs), K(1:dofs, 1:dofs), C(1:dofs, 1:dofs), f, r, 0);
% %Stuck Response
% k_t = 1;%set to 0 to ignore friction
% [ Hhas, Hhtas, phs, Vs, Ds, fdoms ] = LinearFrequencyResponse( M(1:dofs, 1:dofs), K(1:dofs, 1:dofs)+k_t*eye(dofs), C(1:dofs, 1:dofs), f, r, k_t);
% 
% 
% %MASS NORMALISATION
% MR=V'*(M(1:dofs, 1:dofs)*V);%modal mass calculation
% KR=V'*(K(1:dofs, 1:dofs)*V);%modal stiffness calculation
% Vn=V*MR^(0.5);%Normalised eigenvectors
% 
% W_n = sqrt(diag(D)); %pick the maximum frequncy->minimum dt

%% Pattern initialisation
pattern = [1 0 1];
k_t = 1.2;
k_t_v = k_t*pattern';
k_t_v_h = k_t*[pattern 0];
k_t_m = k_t*diag(pattern);

K_t= [ k_t_m -k_t_v];
K_t = [K_t; -k_t_v_h ];
K_t(dofsp, dofsp) = sum(pattern)*k_t;


%% Linear Response with stuck Punch (Dof_punch = dofsp)
K_st = K+K_t;
[V_st,D_st] = eig(K_st,M);%V:eigen vectors, D:eigen values.
d_st=diag(D_st);
W_n_st = sqrt(d_st);


e_st = 4;
r_st = 4;
f_st = zeros(1, dofsp);
f_st(e_st) = 1;


[ Hhast, Hhtast, phst, Vst, Dst, fdomst ] = LinearFrequencyResponse( M,K_st, C, f_st, r_st, k_t);

%% Establish time step parameters
T_n_x = 2*pi/W_n_st(dofsp); % excitation period (max freq min period)
T_n = T_n_x; %pick one for the time steps

dt = T_n/1000;
t_end = 600*T_n; %Number of periods

t = 0:dt:t_end; % total length of time

n = length(t)-1; % number of time steps
t = t(1:n);

%% Determine which special case to use: constant avg. vs. linear accel
method=2; %Const.ave.accel., else linear


if dt/T_n<=0.551 && method ~= 1 % Use linear accel. method - closest to theoretical
    gamma=0.5;
    beta=1/6;
end
if method==1 % Use constant avg. accel. method - unconditionally stable (Example 5.5)
    gamma=0.5;
    beta=0.25;
end

%% Input excitation functions in X (horizontal)
% Input excitation function
P = 1;
unitstep = P*(t>=0);
unitstep2 = P*(t>=1*dt);

%p = t.*unitstep;
p_p=unitstep;
p_p=unitstep-unitstep2;
p_p=p_p(1:n);

w_ext =100; %freq. of excitation
p_p = 1*sin(2*pi*w_ext*t);

P = zeros(dofsp, n);
% P(3, :) = -1*p_p;
% P(e, :) = p_p;
P(dofsp, :) = p_p; %just excite the punch
%P(1, :) = p_p;
%P = repmat(p_x, dofs, 1); %repeat same excitation for all dofs

%% Establish initial conditions @ i=1 

u_x = zeros(dofsp, n);
v_x = zeros(dofsp, n);
a_x = zeros(dofsp, n); 

% u_x(:,1)= 0; % displacement
% v_x(:,1)= 0; % velocity

%%%Calculate Newmark constants

A1_x = M*(1/(beta*dt*dt)) + C*(gamma/(beta*dt));
A2_x = M*(1/(beta*dt))+ C*(gamma/beta-1);
A3_x = M*((0.5/beta)-1) + C*dt*(gamma/(2*beta)-1);

p_hat_x = zeros(dofsp,n);

R_tol = 10^-3;
%% Initialise friction parameters

mu = 0.2;
N0_T = 30; %total normal force
N0 = N0_T/sum(pattern);
%fr_limit = mu*N0;

f_s = zeros(dofsp, n); %restoring force

%% Piezo Pattern and Friction Limit

mu_v = mu*ones(1, dofs);
fr_limit_v = mu_v*N0.*pattern;


%% Run cases
number_case = 1;
switch  number_case
    
    
    
    case 1 %
        close all
        
        %% BASE EXCITATION Horizontal
        t_step = T_n/2;
        P_b = 0;
        u_b_ramp = P_b*t.*(t<=t_step);
        u_b_step = max(u_b_ramp)*(t> t_step);
        u_b = u_b_ramp + u_b_step;
        
        u_b_x = zeros(dofsp, n);
        u_b_x(1, :) = u_b;
        
        v_b_x = zeros(dofsp, n);
        for ii = 1: dofsp
            v_b_x(ii, 1:end-1) = diff(u_b_x(ii, :))./diff(t);
        end
        
        v_b_x(:,n)=v_b_x(:,n-1); %(add one more point at the end)
        
        
        %N = N0*ones(dofs, n); %constant normal load
        
        a_x(:,1) = M\( P(:,1) - C*(v_x(:,1) - v_b_x(:,1)) - K*(u_x(:,1) - u_b_x(:,1)) - f_s(:,1));
 
        for i=1:n-1 %time
            
            p_hat_x(:,i+1) = P(:,i+1)+ A1_x*u_x(:,i)+A2_x*v_x(:,i)+A3_x*a_x(:,i)+ C*v_b_x(:,i+1)+ K*u_b_x(:,i+1);
            
                        
            %% Residual first iteration
            u_x(:,i+1) = u_x(:,i);
            f_s(:,i+1) = f_s(:,i);
                      
            R(:,i+1) = p_hat_x(:,i+1)- K*u_x(:,i+1)- f_s(:,i+1)- A1_x*u_x(:,i+1);
            
            
            K_t= [ k_t_m -k_t_v];
            K_t = [K_t; -k_t_v_h ];
            K_t(dofsp, dofsp) = sum(k_t_v)*k_t;
            
            
            %% Check convergence
            
            while norm(R(:,i+1)) > R_tol
                
                
                K_t= [ k_t_m -k_t_v];
                K_t = [K_t; -k_t_v_h ];
                K_t(dofsp, dofsp) = sum(k_t_v)*k_t;
                
                J = A1_x + K + K_t;
                delta_x = J\R(:,i+1);
                u_x(:,i+1) = u_x(:,i+1) + delta_x;
                
                %Friction in the Membrane Dofs
                f_s(1:dofs,i+1) = f_s(1:dofs,i) + K_t(1:dofs, 1:dofs)*(u_x(1:dofs,i+1)-u_x(1:dofs,i)-u_x(dofsp, i+1)); % within the iteration
                
                               
                %SLIP for each DOF
                
                for ii = 1:dofs
                    fr_limit = fr_limit_v(ii);%for each dofs
                    
                    if k_t ~=0 && abs(f_s(ii,i+1)) > fr_limit    %stick to slip transition
                        %positive
                        if f_s(ii,i+1) > fr_limit
                            
                            ['00a ' num2str(i+1)]
                            
                            f_s(ii,i+1)= fr_limit;
                            
                        end
                        
                        %negative
                        if f_s(ii,i+1) < -fr_limit
                            
                            
                            ['00b ' num2str(i+1)]
                            
                            f_s(ii,i+1)= -fr_limit;
                            
                        end
                        
                    end
                    
                    
                    
                end
                % Friction in Punch is the total of all others
                f_s(dofsp, i+1) = -sum(f_s(1:dofs, i+1));
                
                %update residual
                R(:,i+1) = p_hat_x(:,i+1)- K*u_x(:,i+1)- f_s(:,i+1)- A1_x*u_x(:,i+1);
            end
            
            % udate final values
            v_x(:,i+1) = gamma*(u_x(:,i+1)-u_x(:,i))/(beta*dt) + v_x(:,i)*(1-gamma/beta) + dt*a_x(:,i)*(1-gamma/(2*beta));
            a_x(:,i+1) = (u_x(:,i+1)-u_x(:,i))/(beta*dt*dt) - v_x(:,i)/(beta*dt) -((0.5/beta) -1)*a_x(:,i);
            
        end
        
        %% PLOT
        
        figure (1)
        
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
        
        subplot(2, 3, 4)        
        plot(t, v_x(dofsp,:))
        legend('Punch')
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
        
        subplot(2, 3, 5)        
        plot(t, u_x(dofsp,:))
        legend('Punch')
        xlabel('Time [s]')
        ylabel('Displacement [m]')
        
        subplot(2, 3, 3)
        for i =1:dofs
            txt = ['f_s_', num2str(i)];
            plot(t, f_s(i,:), 'DisplayName', txt)
            hold on
        end
        hold off
        legend show
        xlabel('Time [s]')
        ylabel('Tangential Force [N]')
        
        subplot(2, 3, 6)        
        plot(t, f_s(dofsp,:))
        ylabel('Tangential Force [N]')
        
        
        
        yyaxis right
        hold on
        plot(t, P(dofsp,:))
        legend('Punch Contact force [N]', 'Punch Excitation force [N]')
        xlabel('Time [s]')
        
        
        
        
        figure(5)
        for i =1:dofsp
            subplot(1, dofsp, i)
            txt = ['Hysteresis-DOF_', num2str(i)];
            plot(u_x(i, :), f_s(i, :), 'DisplayName', txt)
            legend show
            xlabel('Displacement [m]')
            ylabel('Tangential Force [N]')
        end
        
        
        figure(8)
        for i =1:dofsp
            txt = ['Residual-DOF_', num2str(i)];
            plot(R(i, :), 'DisplayName', txt)
            hold on
            legend show
        end
        
        
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
        
        figure(10)
        for i=1:dofsp
            
            txt = ['Mode_' , num2str(i)];
            semilogy(fdomst,abs(Hhast(:,i)), 'DisplayName', txt)
            hold on
        end
        grid on
        
        
        semilogy(fdomst,abs(Hhtast), 'DisplayName', ['FRF_e_', num2str(e_st), '_r_', num2str(r_st)])
        legend show
        title('4 DOFs Response - With Punch')
        xlabel('Frequency [Hz]')
          
end
function [ave_curr, conv, dif] = convergence(u, t, n_last, c_c)
%Convergence returns true/1 if the signal converge to the desire parameter
%with moving average
%   x is the signal
%   t is the horizontal axis:time
%   n_last: amount of latest peaks to be considered in the moving
%   average(needs to be a multiple of 3, signal is modulated)
%   T_n: period of the excitation
%   n_per: number of periods considered in the moving average
%   c_c: convergence criterion in percentage

%r = length(u_x(:,1));%number of rows

% a =find(t>t(end)-n_per*T_n);
% 
% t_new = t(a(1):end);
% u_new = u(a(1):end);

t_new = t;
u_new = u;


pks = findpeaks(u_new,t_new); %the values of the peaks
ave_curr = mean(pks(end-n_last+1:end)); %just the last n_last peaks
ave_prev = mean(pks(end-n_last:end-1)); %just the last n_last peaks, 1 before the last peak

dif = 100*((ave_curr-ave_prev)/ave_curr);
conv =0;

if dif <= c_c
    conv =1;
end


%   figure(3)
%     findpeaks(u_x_new(i, :),t_new)
% 
%     text(locs+.02,pks,num2str((1:numel(pks))'))
%     
%     hold on
%     
%     plot(t_new, u_x_new(i,:))
%     
%     plot(t, u_x(i, :))
%     hold off  
    
end


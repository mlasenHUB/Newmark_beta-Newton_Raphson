
number_case = 3;
% case 1: compare with and without plate
% case 2: compare plate for different patterns
% case 3: compare  plate for different patterns (constant k) and without plate
switch  number_case
    case 1
        % compare model with and without plate, with one pattern
        clear all
        close all
        pattern = [1 1 1]; %write here the pattern for the comparison
        %file_with_plate = ['MDF_with_plate_', num2str(pattern), '.mat'];
        file_with_plate = ['MDF_with_plate_k_const_', num2str(pattern), '.mat'];
        file_without_plate = 'MDF_without_plate.mat';
        
        %% Figure 1
        figure(1) %Time response
               
        load(file_with_plate)
        plot(t, u_x, 'LineWidth', 5)
        
        hold on
        %t_end = t(end);
        xlim([0 t(end)])
        
        load(file_without_plate)
        plot(t, u_x, 'LineWidth', 5)
        
        %find the closest time to end
        %[t_end_value,t_end_idx] = min(abs(t-t_end)); %init time index With plate
        legend('u1 with plate', 'u2 with plate','u3 with plate','u4 with plate', 'u5 with plate', 'u4 without plate','u5 without plate'  )
        xlabel('Time [s]')
        
        ylabel('Displacement [m]')
        title(['History ', num2str(pattern),  '-pattern into Plate vs. No Plate'])
        set(gca,'FontSize',25)
        grid on
        
        %% Figure 2
        figure(2)%Hysteresis
        
        load(file_with_plate)
        
        u_mean_plate = mean(u_x(1:dofs, :),1);
        rel_dis = u_x(dp, :)-u_mean_plate;
        plot(rel_dis, f_f(dmm-1, :), 'LineWidth', 3)
        hold on
        plot(u_x(dp, :), f_f(dmm-1, :), 'LineWidth', 3)
        
        load(file_without_plate)
        plot(u_x(1, :), f_f(1, :), 'DisplayName', 'Punch',  'LineWidth', 3)
        
        
        legend('With Plate - Rel.Ave','With Plate - Abs', 'Without Plate' )
        xlabel('Displacement [m]')
        ylabel('Tangential Force [N]')
        title('Comparison responses of the systems with and without plate')
        set(gca,'FontSize',25)
        hold off
        
        %% Figure 3
        figure(3)
        load(file_with_plate)
        t_last = t(end);
        DT = 0.01; %desired time for the last Hysteresis Loops
        t_init = t_last -DT;
        
        %find the closest times in each case and compare the Histeresis Loops
        [init_value_w,t_init_idx_w] = min(abs(t-t_init)); %init time index With plate
        
        %plot
        u_mean_plate = mean(u_x(1:dofs, t_init_idx_w:end),1);
        rel_dis = u_x(dp, t_init_idx_w:end)-u_mean_plate;
        plot(rel_dis, f_f(dmm-1, t_init_idx_w:end), 'LineWidth', 3)
        hold on
        plot(u_x(dp, t_init_idx_w:end), f_f(dmm-1, t_init_idx_w:end), 'LineWidth', 3)
        
        load(file_without_plate)
        [finish_value_wo,t_finish_idx_wo] = min(abs(t-t_last)); % finish time index WithOut plate
        [init_value_wo,t_init_idx_wo] = min(abs(t-t_init)); %init time index With plate
        
        plot(u_x(1, t_init_idx_wo:t_finish_idx_wo), f_f(1, t_init_idx_wo:t_finish_idx_wo), 'DisplayName', 'Punch',  'LineWidth', 3)
        hold off
        
        legend('With Plate - Rel.Ave','With Plate - Abs', 'Without Plate' )
        xlabel('Displacement [m]')
        ylabel('Tangential Force [N]')
        title('Comparison responses of the systems with and without plate')
        set(gca,'FontSize',25)
        hold off
        %% Figure 4
        %State space plots
        load(file_with_plate)
        u_mean_plate = mean(u_x(1:dofs, :),1);
        rel_dis = u_x(dp,:)-u_mean_plate;
        figure(4)
        plot(rel_dis, v_x(dp, :),'DisplayName', 'With Plate', 'LineWidth', 5)
        hold on
        load(file_without_plate)
        plot(u_x(1, :), v_x(1, :),'DisplayName', 'Without Plate', 'LineWidth', 5)
        
        hold off
        legend show
        xlabel('Relative Displacement [m]')
        ylabel('Velocity [m/s]')
        title(['State Space with and without plate; pattern = ', num2str(pattern)])
        set(gca,'FontSize',25)
        
    case 2
        figure (5)
        patterns = [ [0 1 0]; [1 0 1]; [1 1 1]];
        
        for i=1:3
            pattern = patterns(i, :);
            file_with_plate = ['MDF_with_plate_k_const_', num2str(pattern), '.mat'];
            %file_with_plate = ['MDF_with_plate_', num2str(pattern), '.mat'];
            
            load(file_with_plate)
            txtpatterns = ['Pattern: ', num2str(pattern)];
            t_last = t(end);
            DT = 0.01; %desired time for the last Hysteresis Loops
            t_init = t_last -DT;
            
            %find the closest times in each case and compare the Histeresis Loops
            [init_value_w,t_init_idx_w] = min(abs(t-t_init)); %init time index With plate
            %plot
            u_mean_plate = mean(u_x(1:dofs, t_init_idx_w:end),1);
            rel_dis = u_x(dp, t_init_idx_w:end)-u_mean_plate;
            plot(rel_dis, f_f(dmm-1, t_init_idx_w:end), 'DisplayName', txtpatterns,'LineWidth', 5)
            hold on
        end
        legend show
        xlabel('Displacement [m]')
        ylabel('Tangential Force [N]')
        title('Comparison Hysteresis for different patterns')
        set(gca,'FontSize',25)
        hold off
    %% Case 3    
    case 3
        clear all
        close all
        figure (6)
                     
        patterns = [ [0 1 0]; [1 0 1]; [1 1 1]];
        for i=1:3
            pattern = patterns(i, :);
            %file_with_plate = ['MDF_with_plate_k_const_', num2str(pattern), '.mat'];
            file_with_plate = ['MDF_with_plate_', num2str(pattern), '.mat'];
            
            load(file_with_plate)
            txtpatterns = ['Pattern: ', num2str(pattern)];
            t_last = t(end);
            DT = 0.01; %desired time for the last Hysteresis Loops
            t_init = t_last -DT;
            
            %find the closest times in each case and compare the Histeresis Loops
            [init_value_w,t_init_idx_w] = min(abs(t-t_init)); %init time index With plate
            %plot
            u_mean_plate = mean(u_x(1:dofs, t_init_idx_w:end),1);
            rel_dis = u_x(dp, t_init_idx_w:end)-u_mean_plate;
            plot(rel_dis, f_f(dmm-1, t_init_idx_w:end), 'DisplayName', txtpatterns,'LineWidth', 5)
            hold on            
        end
        file_without_plate = 'MDF_without_plate.mat';
        load(file_without_plate)
        [init_value_wo,t_init_idx_wo] = min(abs(t-t_init));
        [finish_value_wo,t_finish_idx_wo] = min(abs(t-t_last)); % finish time index WithOut plate
        
        plot(u_x(1, t_init_idx_wo:t_finish_idx_wo), f_f(1, t_init_idx_wo:t_finish_idx_wo), 'DisplayName', 'No Plate - k_t = 50 N/\mu m','LineWidth', 5)
        
        legend show
        xlabel('Relative Displacement [m] (w.r.t. plate average)')
        ylabel('Tangential Force [N]')
        %title('Comparison Hysteresis for different patterns and No plate -Constant k_t')
        title('Different Patterns on Plate vs. No plate')
        set(gca,'FontSize',25)
        hold off
        
        %% Figure 7
        %state spaces
        
        %State space plots
        figure(7)
        for i=1:3
            pattern = patterns(i, :);
            file_with_plate = ['MDF_with_plate_k_const_', num2str(pattern), '.mat'];
            %file_with_plate = ['MDF_with_plate_', num2str(pattern), '.mat'];
            
            load(file_with_plate)
            txtpatterns = ['Pattern: ', num2str(pattern)];
            u_mean_plate = mean(u_x(1:dofs, :),1);
            rel_dis = u_x(dp,:)-u_mean_plate;
            
            plot(rel_dis, v_x(dp, :),'DisplayName', txtpatterns, 'LineWidth', 5)
            hold on
        
        end
        load(file_without_plate)
        plot(u_x(1, :), v_x(1, :),'DisplayName', 'Without Plate', 'LineWidth', 5)
        
        hold off
        legend show
        xlabel('Relative Displacement [m]')
        ylabel('Velocity [m/s]')
        title('State Space with and without plate for different patterns')
        set(gca,'FontSize',25)
        
end

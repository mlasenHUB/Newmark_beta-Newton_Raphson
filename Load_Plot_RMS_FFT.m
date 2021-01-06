close all
clear all

% %Get filename and path 
%     [fname,pathname] = uigetfile('.csv','Select CSV File to Load, Plot, Compute RMS & FFT');
%     disp([pathname fname])
% 
% %Load CSV
%     tic %start timer
%     data = csvread([pathname fname]); 
%     fprintf('%4.2f seconds - Time to Load Data\n',toc)
load('E:\From Time Integration\MDF_without_plate_60N_tang_Load.mat')
    
%Determine variables and Display size
    %[N,m] = size(data);
    t=t(4000:end);
    N= length(t);
    x=u_x(1, 4000:end);
%     N= 2000;
%     t = t(end-N+1:end);
%     x= u_x(1, end-N+1:end);
    Fs = 1/(t(2)-t(1));
    fprintf('%12.0f data points\n',N)

%Plot Data
    tic %start timer
    figure(1)
    plot(t,x)
    xlabel('Time (s)');
    ylabel('Accel (g)');
%    title(fname);
    grid on;
    fprintf('%4.2f seconds - Time to Plot Data\n',toc)
    
%Determine RMS and Plot
    tic %start timer
    %% figure 1
    w = floor(Fs); %width of the window for computing RMS
    steps = floor(N/w); %Number of steps for RMS
    t_RMS = zeros(steps,1); %Create array for RMS time values
    x_RMS = zeros(steps,1); %Create array for RMS values
    for i=1:steps
        range = ((i-1)*w+1):(i*w);
        t_RMS(i) = mean(t(range));
        x_RMS(i) = sqrt(mean(x(range).^2));  
    end
    %% figure 2
    figure(2)
    plot(t_RMS,x_RMS)
    xlabel('Time (s)');
    ylabel('RMS Accel (g)');
    %title(['RMS - ' fname]);
    grid on;
    fprintf('%4.2f seconds - Time to Compute RMS and Plot\n',toc)    
    %% figure 3
%Determine FFT and Plot
    tic 
    freq = 0:Fs/length(x):Fs/2; %frequency array for FFT
    xdft = fft(x); %Compute FFT
    xdft = 1/length(x).*xdft; %Normalize
    xdft(2:end-1) = 2*xdft(2:end-1);
    
    figure(3)
    semilogy(freq,abs(xdft(1:floor(N/2)+1)), 'DisplayName', 'No Plate', 'LineWidth', 5)
    xlabel('Frequency (Hz)');
    ylabel('Amplitude (m)');
%    title(['FFT - ' fname]);
    grid on;
    hold on 
    fprintf('%4.2f seconds - Time to Compute FFT and Plot\n',toc)
    
    %% load 010
    load('E:\From Time Integration\test_0  1  0.mat', 't', 'u_x', 'f_f', 'W_n_st')
    t=t(4000:end);
    N= length(t);
    x=u_x(4, 4000:end);
    Fs = 1/(t(2)-t(1));
    
    freq = 0:Fs/length(x):Fs/2; %frequency array for FFT
    xdft = fft(x); %Compute FFT
    xdft = 1/length(x).*xdft; %Normalize
    xdft(2:end-1) = 2*xdft(2:end-1);
    
    figure(3)
    semilogy(freq,abs(xdft(1:floor(N/2)+1)), 'DisplayName', '010', 'LineWidth', 5)
    
    %% 101
    %% load 010
    load('E:\From Time Integration\test_1  0  1.mat', 't', 'u_x', 'f_f', 'W_n_st')
    t=t(4000:end);
    N= length(t);
    x=u_x(4, 4000:end);
    Fs = 1/(t(2)-t(1));
    
    freq = 0:Fs/length(x):Fs/2; %frequency array for FFT
    xdft = fft(x); %Compute FFT
    xdft = 1/length(x).*xdft; %Normalize
    xdft(2:end-1) = 2*xdft(2:end-1);
    
    figure(3)
    semilogy(freq,abs(xdft(1:floor(N/2)+1)), 'DisplayName', '101', 'LineWidth', 5)
    %% 111
    %% load 010
    load('E:\From Time Integration\test_1  1  1.mat', 't', 'u_x', 'f_f', 'W_n_st')
    t=t(4000:end);
    N= length(t);
    x=u_x(4, 4000:end);
    Fs = 1/(t(2)-t(1));
    
    freq = 0:Fs/length(x):Fs/2; %frequency array for FFT
    xdft = fft(x); %Compute FFT
    xdft = 1/length(x).*xdft; %Normalize
    xdft(2:end-1) = 2*xdft(2:end-1);
    
    figure(3)
    semilogy(freq,abs(xdft(1:floor(N/2)+1)), 'DisplayName', '111', 'LineWidth', 5)
    
    xlim([0 3000])
    set(gca,'FontSize',40)
    legend show
    hold off
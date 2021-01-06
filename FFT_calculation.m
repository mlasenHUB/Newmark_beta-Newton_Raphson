%% Calculate the FFT  of the Lumped Model results
close all
%load('E:\From Time Integration\MDF_with_plate_extracted_kt_from_simulations_C6_and_diff_k_plate0  1  0.mat', 't', 'u_x', 'f_f', 'W_n_st')
%load('E:\From Time Integration\test_0  1  0.mat', 't', 'u_x', 'f_f')
load('E:\From Time Integration\MDF_without_plate_60N_tang_Load.mat')

y = u_x(1, :);
L = length(y);
NFFT = 2^nextpow2(L); % Next power of 2 from length of y
dt = t(2) -t(1);
Fs = 1/dt; %sampling frequency
Y = fft(y,NFFT)/L;
f = Fs/2*linspace(0,1,NFFT/2+1);%frequency vector




figure
plotLength=round(length(f)/2);
semilogy(f(1:plotLength),2*abs(Y(1:plotLength))) 

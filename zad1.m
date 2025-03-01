close all;
clear;clc;
%% signal parameters
f_sin = 1e9;%frequency of modulated signal
c = physconst('LightSpeed');
lambda = c/f_sin; %wavelenght
pulseDuration = 10e-6;%time of signal transmission
B = 2e6; %bandwidth            
mi=B/pulseDuration; %modulation coefficient
PRI = 1/100;%frequency of azimuth sampling
PRF=1/PRI;
TotalTime = 0.5;         % total time of meas , unused
f_sample = 60e6;        % range sample frequency

% param calculation 
% number of pulses (azimuth samples)
N_pulses = floor(TotalTime/PRI);  
%sample time 
t_sample = 1/f_sample;
%range resolution
dR = c/(2*f_sample);

%% geometry params
%параметры спутника - расстояние , скорость 
satelite_location = [0, 0, 100e3];  
V_satelite = [0, 150, 0];     
%параметры цели 
object_location = [1000,5000, 0];   
%max distance from radar to target 
Rmax=120e3; % previously it was 30e3 , but its not enough if satellite altitude is 100km

%closest range
R0_nonorm=[object_location(1)-satelite_location(1),0,object_location(3)-satelite_location(3)];
R0=norm(R0_nonorm);
%num of range samples
N_range = floor(Rmax/dR);
% matrix for initial raw data N_pulses columns and N_range rows
%матрица для данных - Колво импульсов по вертикали(колонки) и расстояние по
%горизонтали (ряды)
raw_data = zeros(N_range, N_pulses);

%% calculating transmitted signal
t_chirp = 0:t_sample:(pulseDuration-t_sample);
chirp = exp(2i*pi*(mi*t_chirp-mi*pulseDuration).*t_chirp/2);

%% precalculation 
% used to calculate sattelite path timings
pulse_times = (0:N_pulses-1) * PRI;
% calculating whole path of satellite
satellite_positions = satelite_location + V_satelite .* pulse_times';
satellite_positions = reshape(satellite_positions, N_pulses, 3);

f_doppler=zeros(1,N_pulses);
%%calculation per every azimuth sample
for i = 1:N_pulses
    % distance calculation
    distance = object_location - satellite_positions(i,:);
    d_norm = norm(distance);
    % doppler shift calculation
    V_rel = dot(V_satelite, distance/d_norm);
    %azimuth chirp rate
    
    %doppler frequency
    f_doppler(i) = (2*V_rel*f_sin)/c;
    %calculation range cell number 
    range_idx = floor(d_norm/dR);
    % calculating phase 
    phase = -4*pi*d_norm/lambda - 2*pi*f_doppler(i)*pulse_times(i);
    raw_data(range_idx, i) = exp(1i * phase);
end

%% raw data
figure;
imagesc(abs(raw_data));
xlabel('Azimuth freq');
ylabel('Range samples');
title('Raw SAR data');
colorbar;

mathed_filter=conj(flipud(chirp));%flip just to turn horizontal to vertical

%compression first approach
F=fftshift(fft(chirp,N_range));
window=hamming(N_range);
F=F'.*window;
R_cmp=zeros(size(raw_data));
for i=1:N_pulses
    R=fft(raw_data(:,i),N_range);
    R_cmp(:,i)=R.*F;
    R_cmp(:,i)=ifft(R.*F);
end

%% azimuth fft
raw_az_fft=fftshift(fft(R_cmp',[],1),1);
figure;
imagesc(abs(raw_az_fft));
title('range doppler map');
xlabel('Range samples');
ylabel('Azimuth freq');



%% rcmc
delta_R = (lambda^2 * f_doppler.^2 * R0) / (8 * V_satelite(2)^2);
delta_samples = delta_R / dR;  
raw_RCMC= zeros(size(raw_az_fft));
%{
%first approach from book
[Frg_2D, dR2D] = meshgrid(Frg,dR);
G = exp(1j * 4 * pi * Frg_2D .* dR2D / c);
raw_RCMC=raw_az_fft.*G;
%}
%another approach from book
for i= 1:N_pulses
  shift = delta_samples(i);
    
  % calcof bins to shift
  range_bins = 1:N_range;
  shifted_bins = range_bins - shift;
  
  % interpolation for shifting range samples
  raw_RCMC(i,:) = interp1(range_bins, raw_az_fft(i,:), shifted_bins, 'linear', 0);
end
figure;
imagesc(abs(raw_RCMC));
title('rcmc');
xlabel('Range samples');
ylabel('Azimuth freq');

%% azimuth compression 
Kaz=(2*V_satelite(2)^2)/(lambda*R0);
t_az=(-N_pulses/2:(N_pulses-1)/2)/PRF;
h_az=exp(-1i*pi*Kaz*t_az.^2);
window_azimuth = hamming(N_pulses)';
h_az = h_az .* window_azimuth;
for i = 1:N_range
    raw_RCMC(:,i) = raw_RCMC(:,i) .* h_az';
end
%% azimuth ifft

sar_image=ifft(raw_RCMC,[],2);
figure;
imagesc(abs(sar_image));
xlabel('Azimuth freq');
ylabel('Range samples');
title('SAR data after rda');
colorbar;




close all;
clear;clc;

%% signal parameters
f_sin=1e9;%1GHz
c=physconst('LightSpeed');
lambda=c/f_sin;
A=1;%1V unused
pulseDuration=10e-6; %100us
B=2e6;%bandwidth in hz
%modulation index
mi=B/pulseDuration;
%max range
Rmax=100e3;
%pulse range time
PRI=0.1e-3;%1 ms
%total time of simulation
TotalTime=10e-3;%10ms
%sample time & frequency
f_sample=4e6;%50MHz
t_sample=1/f_sample;
%range resolution
dR=c/(2*f_sample);
%% point parameter - location & RCS
object_location=[100,100,0]; %x,y,z in km
object_RCS=0.5;
%satelite parameters - location & movement speed
satelite_location=[0,0,100e3]; %x,y,z in km
V_satelite=[0,1000,0]; %1km/s
%satellite is moving in the direction of y axis 

%% signal generation
%time vector for whole time, unused for now
t=0:t_sample:TotalTime;
%time vector for a single pulse
timp=0:t_sample:(pulseDuration-t_sample);

%signal 
%t_signal=0:1/f_sample:pulseDuration;
%signal=A*sin(2*pi*f_sin*t_signal);
%signal is then used in a loop to simulate echo
signal=exp(2i*pi*(mi*timp-mi*pulseDuration).*timp/2);
%imp_time=0:PRI:TotalTime;

%% calculation of values needed to future loop
% number of samples
N=TotalTime/PRI;
%calculating num of samples per period
N_sample=floor(Rmax/dR);
%to save the number of samples in one pulse
sample_per_pulse=N_sample; 
%scaling to number of measurements
N_sample=N_sample*N;
tabdelt=zeros(1,N_sample);
%satellites travaled distance for 10ms - 10km
satellite_path=TotalTime*V_satelite;
%% echo , freq and phase calculation
f_doppler=zeros(1,N);
s=zeros(1,N);
for i=1:(N)
    %distance between satellite and object
    distance=object_location-(satelite_location+V_satelite*i*PRI);
    %distance between satellite and object in i-th moment
    d_norm=norm(distance);
    d_i=distance/d_norm;
    %recalculation of sample number for each norm distance
    nrprob=floor(d_norm/dR);
    %relative velocity in i-th moment
    V_rel=V_satelite*d_i';
    %doppler shift in i-th moment , needed to 
    f_doppler(i)=(2*V_rel*f_sin)/c;
    %echo recalculation 
    tabdelt(nrprob+(i-1)*sample_per_pulse)=tabdelt(nrprob+(i-1)*sample_per_pulse) + object_RCS*exp( -4i*pi*d_norm/lambda  +  2i*pi*f_doppler(i) );
    temp=conv(tabdelt,signal);
    sygnOdb=temp(1:N_sample);

end
%% compression
sygnOdb_matrix=reshape(sygnOdb,sample_per_pulse,N)';
range_compressed=zeros(N,sample_per_pulse);
for i=1:N
    range_compressed(i,:)=conv(sygnOdb_matrix(i,:),signal,'same');
end

% tu nic waznego , tylko wykresy
figure;
imagesc(abs(range_compressed));
range_axis = (0:sample_per_pulse-1) * dR; % in meters
azimuth_axis = (0:N-1) * PRI * 1000; % Approximate azimuth in meters
figure;
imagesc(range_axis, azimuth_axis, abs(range_compressed));

%filtracja
window=hamming(N);
range_compressed=range_compressed.*window;

%azimuth fft
range_doppler_map=fftshift(fft(range_compressed,[],1),1);

%azimuth ifft
sar_image=ifft(range_doppler_map,[],1);
figure;
imagesc(abs(sar_image));
%% plotting
%plot
%figure;
%plot(f_doppler);
%this is a doppler shift for each measurement . Smallest one on the 50 km , where the target is
%figure;
%plot(real(s));
%hold on;
%plot(imag(s));
%figure;
%plot(abs(sygnOdb),'b');
%figure;
%plot(real(sygnOdb),'b');
%hold on;
%plot(imag(sygnOdb),'r');

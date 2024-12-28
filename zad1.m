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
Rmax=30e3;
%pulse range time
PRI=0.1e-3;%1 ms
%total time of simulation
TotalTime=100e-3;%10ms
%sample time & frequency
f_sample=4e6;%50MHz
t_sample=1/f_sample;
%range resolution
dR=c/(2*f_sample);
%% point parameter - location & RCS
object_location=[20,20,0]; %x,y,z in km
object_RCS=0.5;
%satelite parameters - location & movement speed
satelite_location=[0,0,100]; %x,y,z in km
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
ileprob=floor(Rmax/dR);
%to save the number of samples in one pulse
probe_per_pulse=ileprob; 
%scaling to number of measurements
ileprob=ileprob*N;
tabdelt=zeros(1,ileprob);
%satellites travaled distance for 10ms - 10km
satellite_path=TotalTime*V_satelite;
%% echo , freq and phase calculation
f_doppler=zeros(1,N);
s=zeros(1,N);
for i=1:(N-1)
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
    tabdelt(nrprob+i*probe_per_pulse)=tabdelt(nrprob+i*probe_per_pulse) + object_RCS*exp( -4i*pi*d_norm/lambda  +  2i*pi*f_doppler(i) );
    temp=conv(tabdelt,signal);
    sygnOdb=temp(1:ileprob);
    %s(i)=exp((-4i*pi*d_norm)/lambda);
end

%% plotting
%plot
figure;
plot(f_doppler);
%this is a doppler shift for each measurement . Smallest one on the 50 km , where the target is
%figure;
%plot(real(s));
%hold on;
%plot(imag(s));
figure;
plot(abs(sygnOdb),'b');
figure;
plot(real(sygnOdb),'b');
hold on;
plot(imag(sygnOdb),'r');


clear;clc;

%% signal parameters
f_sin=1e6;%1MHz
lambda=physconst('LightSpeed')/f_sin;
A=1;%1V
pulseDuration=100e-6; %100us
PRF=1e-3;%1 ms
TotalTime=100e-3;%10ms
f_sample=50e6;%50MHz

%% point parameter
object_location=[50,50,0]; %x,y,z in km
satelite_location=[0,0,100]; %x,y,z in km
V_satelite=[0,1000,0]; %1km/s
%satellite is moving in the direction of y axis 

%% signal generation
%time vector
t=0:1/f_sample:TotalTime;

%signal 
t_signal=0:1/f_sample:pulseDuration;
signal=A*sin(2*pi*f_sin*t_signal);

imp_time=0:PRF:TotalTime;

pulse=pulstran(t,imp_time,signal,f_sample);

%%
%plot
plot(t,pulse);

%%
% number of samples
N=TotalTime/PRF;
% N samples of reflected signal

%satellites travaled distance for 10ms - 10km
satellite_path=TotalTime*V_satelite;
%line of sight 
%LOS=[satelite_location(1)-object_location(1),satelite_location(2)-object_location(2),satelite_location(3)-object_location(3)];
%LOS_1=LOS/norm(LOS);
%relative velocity along the line of sight
%V_rel=V_satelite*LOS_1';

f_doppler=zeros(1,N);
for i=1:N
    %distance between satellite and object
    d=satelite_location-object_location;
    %distance between satellite and object in i-th moment
    d_norm=norm(satelite_location+V_satelite*i*PRF-object_location);
    d_i=d/d_norm;
    %relative velocity in i-th moment
    V_rel=V_satelite*d_i';
    %V_rel_i=V_satelite*i*PRF*LOS_1';
    %doppler shift in i-th moment
    %f_doppler(i)=f_sin*(d_i-d)/d-V_rel_i/physconst('LightSpeed');
    f_doppler(i)=(2*V_rel*f_sin)/physconst('LightSpeed');

end

%plot
figure;
plot(f_doppler);
%this is a doppler shift for each measurement . Smallest one on the 50 km , where the target is

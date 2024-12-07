
clear;clc;

%% signal parameters
f_sin=1e9;%1MHz
c=physconst('LightSpeed');
lambda=c/f_sin;
A=1;%1V
pulseDuration=100e-6; %100us
PRI=0.1e-3;%1 ms
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

imp_time=0:PRI:TotalTime;

pulse=pulstran(t,imp_time,signal,f_sample);

%%
%plot
plot(t,pulse);

%%
% number of samples
N=TotalTime/PRI;
% N samples of reflected signal

%satellites travaled distance for 10ms - 10km
satellite_path=TotalTime*V_satelite;
%line of sight 
%LOS=[satelite_location(1)-object_location(1),satelite_location(2)-object_location(2),satelite_location(3)-object_location(3)];
%LOS_1=LOS/norm(LOS);
%relative velocity along the line of sight
%V_rel=V_satelite*LOS_1';

f_doppler=zeros(1,N);
s=zeros(1,N);
for i=1:N
    %distance between satellite and object
    distance=object_location-satelite_location+V_satelite*i*PRI;
    %d=satelite_location-object_location;
    %distance between satellite and object in i-th moment
    d_norm=norm(distance);
    d_i=distance/d_norm;
    %relative velocity in i-th moment
    V_rel=V_satelite*d_i';
    
    %doppler shift in i-th moment
    f_doppler(i)=(2*V_rel*f_sin)/c;
    s(i)=exp((1i*4*pi*d_norm)/lambda);
end

%plot
figure;
plot(f_doppler);
%this is a doppler shift for each measurement . Smallest one on the 50 km , where the target is
figure;
plot(real(s));
hold on;
plot(imag(s));

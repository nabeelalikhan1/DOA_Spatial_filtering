%% computing MSE curve vs. SNR in under-determined scenario
%In this code we assume that the user add in this path the TFSAP toolboox.
%Please cite the following reference:

close all;
clear all;
addpath('D:\tfsa_5-5\windows\win64_bin');
%addpath('C:\Users\DR. Nabeel\Desktop\IF estimation using connectivity\DOA-estimation-of-intersecting-components-master\DOA-estimation-of-intersecting-components-master');

% NUmber of simulation runs

% Signal model%
%addpath 'D:\tfsa_5-5\windows\win64_bin\'
index=0;
n=0:127;
% Number of SOurces
% NUmber of components
s2=exp(2*pi*1i*(1*0.1*n+1*0.3*n.^3/(128*128*3)));

s1=exp(2*pi*1i*(0.5*n-1*.3*n.^3/(128*128*3)));

s3=exp(2*pi*1i*(0.4*n-0*0.3*n.^3/(128*128*3)));
s5=exp(2*pi*1i*(0.435*n-0.3*n.^3/(128*128*3)));
s6=exp(2*pi*1i*(0.05*n+0.25*n.^3/(128*128*3)));
s7=exp(2*pi*1i*(0.1*n+0.15*n.^3/(128*128*3)));
SampFreq=128;
t = 0:1/SampFreq:1-1/SampFreq;
s8=exp(1*1i*(2*pi*(SampFreq*t/4 +2*sin(4*pi*t))));
s9=exp(1*1i*(2*pi*(SampFreq*t/4 -2*sin(4*pi*t))));
s4=exp(2*pi*1i*(0.15*n-0.0*n.^3/(128*128*3)));
s = [(s1.') (s2.') (s8.') (s9.') ];
%s = [(s8.')  (s9.')];
s = [(s1.') (s2.') (s5.') (s6.') ];

%s = [(s1.') (s2.') ];

% set mixing matrix A
theta = [10,25,40,60]*pi/180;   % sensor separation angles in radians
%theta = [0,10,20]*pi/180;   % sensor separation angles in radians

%theta = [15,25]*pi/180;   % sensor separation angles in radians
n_sources=length(theta);
N_sensors=8;

M=N_sensors;% Number of sensors

% Signal generation
s_orig=s;
% set mixing matrix A
i_snr=0;
A = exp(1j*pi*[0:M-1].'*sin(theta));  % mixing matrix A
    
X = A*s.';                             % mixed source
theta9=round(theta *180/pi);
% genera/2te noise
SNR=-10;
sigma = 10^(-SNR/20);
w = sigma*(randn(M,length(n)) + 1j*(randn(M,length(n))))/sqrt(2); % noise
Xo=X+w;
DOA=STFD_SPATIAL_FILTERIN_Illustration(Xo,n_sources,N_sensors);


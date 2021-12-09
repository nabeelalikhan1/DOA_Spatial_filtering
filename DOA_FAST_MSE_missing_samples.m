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
for NSS=16:16:80
SNR=0;
for sim_c=1:100
    
%X = A*s.';   
B=rand(size(A));
X=(1.*A)*s.';
% mixed source
theta9=round(theta *180/pi);
% genera/2te noise
sigma = 10^(-SNR/20);
w = sigma*(randn(M,length(n)) + 1j*(randn(M,length(n))))/sqrt(2); % noise
%X=X+w;
for i=1:M
p=randperm(128);

    X(i,p(1:NSS))=0;
end
Xo=X;

for i=1:n_sources
    
P = TMMUSIC(X*X', 2, N_sensors, n_sources-i+1, 1, (0:1:90)');
    [V,ind]=max(P);
    %figure;plot(P)
   DOA_M(i)=ind-1;
        AA = exp(-1j*pi*[0:M-1].'*sin(round(DOA_M(i))*pi/180)).';
    x = (AA)*Xo/N_sensors;
 X=X-(AA')*x;%*abs(sum(AA*X)).^2;
end




[IFF,ss] = Multi_Sensor_FASTEST_IF(Xo,N_sensors,65, n_sources, 3,50,0,0,1,length(X));

for iii=1:n_sources
    for jjj=1:N_sensors
        a(jjj,:)=ss(jjj,iii,:);
    end
    theta1=0:1:180;
    
    p=TMMUSIC(a*a', 2, N_sensors, 1, 1, (0:1:90)');
    % figure; plot(p);
    [x,y]=max(p);
    y1(iii)=y(1);
end
%[ss,IFF] = multi_sensor_source_separation_spatial_diversity(Xo, 2,N_sensors,0);
ss= multi_sensor_source_separation(Xo, n_sources, 3,N_sensors);

[aa,bb,~]=size(ss);
clear y2;
for iii=1:bb%n_sources
    for jjj=1:aa%N_sensors
        a(jjj,:)=ss(jjj,iii,:);
    end
    theta1=0:1:90;
    
    p=TMMUSIC(cov(a'), 2, N_sensors, 1, 1, theta1');
    [x,y]=max(p);
    y2(iii)=y(1);
    
end
if length(y2)>length(theta)
            y2=y2(1:4);
        elseif length(y2)<length(theta)
            y2(length(y2):length(theta))=0;
        end
% [IFF,ss] = Multi_Sensor_FASTEST_IF_spatial(Xo,N_sensors,65, n_sources, 3,50,0,0,1,length(X));
%     theta1=0:1:180;
% 
% for iii=1:n_sources
%     for jjj=1:N_sensors
%         a(jjj,:)=ss(jjj,iii,:);
%     end
%     
%     p=TMMUSIC(a*a', 2, N_sensors, 1, 1, (0:1:90)');
%     % figure; plot(p);
%     [x,y]=max(p);
%     y2(iii)=y(1);
% end
%y1=y1-90;
%DOA=FAST_DOA_Spatial_filtering(Xo,N_sensors,75, n_sources, 100,length(X));
DOA=STFD_SPATIAL_FILTERIN(Xo,n_sources,N_sensors);

DOA_M=sort(DOA_M);
y1=sort(y1);
y2=sort(y2);
%y2=sort(y2);
DOA=sort(DOA);
theta11=sort(theta*180/pi);
y1=y1/max(theta11);
y2=y2/max(theta11);

%y2=y2/max(theta11);
DOA=DOA/max(theta11);
DOA_M=DOA_M/max(theta11);
theta11=theta11/max(theta11);

mse_MUSIC(sim_c)=mean(abs(DOA_M-theta11));
mse_spatial_Filtering(sim_c)=mean(abs(DOA-theta11));
mse_FAST_DOA(sim_c)=mean(abs(y1-theta11));
mse_connect(sim_c)=mean(abs(y2-theta11));

%mse_FAST_DOA_coherent(sim_c)=mean(abs(y2-theta11));

end
i_snr=i_snr+1;
mse_MUSIC1(i_snr)=mean(mse_MUSIC)
mse_spatial_Filtering1(i_snr)=mean(mse_spatial_Filtering)
mse_FAST_DOA1(i_snr)=mean(mse_FAST_DOA)
mse_c(i_snr)=mean(mse_connect)
%mse_FAST_DOA_coherent1(i_snr)=mean(mse_FAST_DOA_coherent);

end

mse_MUSIC1
mse_spatial_Filtering1
mse_FAST_DOA1
%mse_FAST_DOA_coherent1

%
% [ss,IFF] = multi_sensor_source_separation_spatial_diversity(X, 3,N_sensors,0);
% [aa,bb,~]=size(ss);
% clear y1;
% clear a;
% theta1=-90:90;
% for iii=1:bb%n_sources
%     for jjj=1:aa%N_sensors
%         a(jjj,:)=ss(jjj,iii,:);
%     end
%
%     p=TMMUSIC(cov(a'), 2, N_sensors, 1, 1, theta1');
%     [x,y]=max(p);
%     y1(iii)=y-91;
%
% end
%
%

plot(16:16:80,20*log10(mse_MUSIC1),'r','LineWidth',3);
hold on;plot(16:16:80,20*log10(mse_spatial_Filtering1),'b:','LineWidth',3);
hold on;plot(16:16:80,20*log10(mse_FAST_DOA1),'k-','LineWidth',3);
hold on;plot(16:16:80,20*log10(mse_c),'g-','LineWidth',3);

xlabel('Number of Missing Samples')
ylabel('Mean Absolute Error (dB)')
legend('Conventaional time-domain MUSIC algorithm with spatial filtering','The Proposed Method','DOA estimation based on IF Estimation','Connectivity based Method'); 
axis([16 80 -40 -20])



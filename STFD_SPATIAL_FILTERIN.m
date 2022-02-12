%% computing MSE curve vs. SNR in under-determined scenario
%In this code we assume that the user add in this path the TFSAP toolboox.
%Please cite the following reference:
function DOA_P=STFD_SPATIAL_FILTERIN(X,n_sources,N_sensors)

Xo=X;
%%%%%%%%  BSS code  ends
for i=1:n_sources
   %  D = mtfd(X, 'CKD',0.15,0.15);
  % D_avg = zeros(length(X), length(X));
  % for mm = 1:N_sensors, D_avg = D{mm,mm} + D_avg; end
  % D_avg=real(D_avg);
    
 %[D_avg,D]=RSTFD(D,N_sensors);
   
    [D_avg,D,~]= SADTFD_new(X,3,8,64);
        %[D_avg,D,~]= SADTFD_new(X,2,20,64);

    thr=0.1*max(D_avg(:))*6;
    thr=0.1*max(D_avg(:))*9;
    
    %peaks = peak_tfd(tfd, thr);
    peaks=zeros(size(D_avg));
    D_avg(D_avg<0)=0;
    peaks(D_avg>=thr)=1;
    D_avg=D_avg.*peaks;

    orient=zeros(size(D_avg))-1000;
    for ii=1:length(D_avg)
        for jj=1:length(D_avg)
            if peaks(ii,jj)==1
                Ds=zeros(N_sensors,N_sensors);
                for mm = 1:N_sensors
                    for nn=1:N_sensors
                        Ds(mm,nn)=D{mm,nn}(ii,jj);
                    end
                end
                [~,d]=eig(Ds);
                d=sort(abs(diag(d)));

                P = TMMUSIC(Ds, 2, N_sensors, 1, 1, (0:1:90)');
                [~,orient(ii,jj)]=max(P);
                if (abs(d(end-1)))>0.5*abs(d(end))
                    %P = tf_music(Ds, 1, N_sensors, 2, 1, (-90:1:90)');
                %      D_avg(ii,jj)=0;
                end
            end
        end
    end
    Thr2=4;
     % figure; imagesc(D_avg);
   DOA_mean= orient(D_avg==max(D_avg(:)));
%    orient=orient-91;
   [ii,jj]=find(and(orient<=DOA_mean+Thr2,orient>=DOA_mean-Thr2));
                Ds=zeros(N_sensors,N_sensors);
    for k=1:length(ii)
        for mm = 1:N_sensors
            for nn=1:N_sensors
                Ds(mm,nn)=Ds(mm,nn)+D{mm,nn}(ii(k),jj(k));
            end
        end
    end
    [~,d]=eig(Ds);
    %P = tf_music(Ds, 1, N_sensors, 2, 1, (0:1:180)');
    P = TMMUSIC(Ds, 2, N_sensors, 1, 1, (0:1:90)');
    
    [v,ind]=  max(P);
    DOA_mean=ind-1;

    DOA_P(i)=DOA_mean;
    % Code for spatial filtering
    AA = exp(-1j*pi*[0:N_sensors-1].'*sin(round(DOA_mean)*pi/180)).';
    x = (AA)*Xo/N_sensors;
    
   X=X-(AA')*x;%*abs(sum(AA*X)).^2;
end
end
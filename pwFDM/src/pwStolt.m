function [BDcomp,BDall] = pwStolt(RFdata,settings)

%--------------------------------------------------------------------------
%-- Modified Stolt's method: RF data, PW centered at x = 0
%-- (Note: Beamformed data is compounded in Kz-x domain)
%-- Authors: M.Albulayli and D.Rakhmatov (daler@ece.uvic.ca)
%-- Institution: ECE Department, University of Victoria
%-- Date: 2018/04/11
%--------------------------------------------------------------------------


[Nt,Nx,Nf] = size(RFdata);
Nz = Nt;


%% Input specifications

c = settings.c0;            %-- speed of sound (m/s)
fs = settings.fs;           %-- sampling frequency (Hz)
theta = settings.PWangle;   %-- PW emission angles (rad)

t = settings.t;             %-- time: t-axis (s)
z = settings.z;             %-- imaging depth: z-axis (m)
x = settings.x;             %-- lateral coordinates: x-axis (m)

dz = mean(diff(z));         %-- depth spacing (m)
dx = mean(diff(x));         %-- lateral spacing (m)

v = c/2;                    %-- ERM velocity assumption (m/s)


%% Configure 2-D FFT grid

NxFFT = settings.NxFFT;     %-- FFT points for x-axis
NtFFT = settings.NtFFT;     %-- FFT points for t-axis

Kx0 = (-NxFFT/2:NxFFT/2-1)/dx/NxFFT;    %-- spatial frequencies
F0 = (0:NtFFT/2-1)*fs/NtFFT;            %-- positive temporal frequencies

[Kx,F] = meshgrid(Kx0,F0);


%% Initial (zero-valued) beamformed data in Kz-x and z-x domains

D_Kz_x_comp = zeros(NtFFT/2,Nx);    %-- compounded data in Kz-x domain

BDcomp = zeros(Nz,Nx);              %-- compounded data in z-x domain
BDall = zeros(Nz,Nx,Nf);            %-- all individual frames


%% Stolt's mapping for each plane wave:

for j = 1:Nf

    D_t_x = RFdata(:,:,j);
    
    D_F_x = fft(D_t_x,NtFFT);   
    D_F_x = D_F_x(1:NtFFT/2,:);         %-- extract positive-F spectrum    
    D_F_Kx = fft(D_F_x,NxFFT,2);
    D_F_Kx = fftshift(D_F_Kx,2);        %-- shift zero-Kx to center
    
    %-- remove evanescent parts in F-Kx domain
    evanescent = (F.^2 <= c^2*Kx.^2);
    D_F_Kx(evanescent) = 0;    
    
    %-- determine frequency interpolation points
    Kz = F/v;                           %-- Kz corresponds to z = t*c/2
    Knew = Kz.*(1 + (Kx./Kz).^2);
    Knew(isnan(Knew)) = 0;
    compressZ = (1+cos(theta(j)))/2;
    Fnew = v*Knew/compressZ;
    
    %-- map from F-Kx to Kz-Kx domain
    D_Kz_Kx = zeros(NtFFT/2,NxFFT);
    for i = 1:NxFFT
        D_Kz_Kx(:,i) = interp1(F(:,i),D_F_Kx(:,i),Fnew(:,i),'linear',0);
    end
    
    %-- scale in Kz-Kx domain
    negative = (Kz.^2 <= Kx.^2);
    scaler = v*(1 - (Kx./Kz).^2)/compressZ;
    scaler(negative) = 0;
    D_Kz_Kx = scaler.*D_Kz_Kx;
    
    %-- transform back to Kz-x domain
    D_Kz_Kx = ifftshift(D_Kz_Kx,2);       
    D_Kz_x = ifft(D_Kz_Kx,NxFFT,2);    
    D_Kz_x = D_Kz_x(:,1:Nx);
    
    %-- compensate for plane-wave angle in Kz-x domain
	if theta(j)~=0       
        for i = 1:Nx
            zshift = x(i)*tan(theta(j))/2;
            D_Kz_x(:,i) = D_Kz_x(:,i).*exp(1i*2*pi*zshift*Kz(:,i));
        end
    end
    
    %-- compound data in Kz-x domain (positive spectrum only)
    D_Kz_x_comp = D_Kz_x_comp + D_Kz_x;
    
    %-- transform current frame into z-x domain (optional record keeping)
    D_z_x = ifft(D_Kz_x,NtFFT,'symmetric');     %-- symmetric negative spectrum   
    D_z_x = D_z_x(1:Nz,:);

    BDall(:,:,j) = D_z_x;
    
end

%-- transform compounded data into z-x domain
BDcomp = ifft(D_Kz_x_comp,NtFFT,'symmetric');	%-- symmetric negative spectrum
BDcomp = BDcomp(1:Nz,:);

end

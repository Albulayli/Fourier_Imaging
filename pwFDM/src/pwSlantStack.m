function [BDcomp,BDall] = pwSlantStack(RFdata,settings)

%--------------------------------------------------------------------------
%-- Modified slant-stack migration: RF data, PW centered at x = 0
%-- (Note: Beamformed data is compounded in z-x domain)
%-- Authors: M.Albulayli and D.Rakhmatov (daler@ece.uvic.ca)
%-- Institution: ECE Department, University of Victoria
%-- Date: 2018/05/12
%--------------------------------------------------------------------------


[Nt,Nx,Nf] = size(RFdata);
Nz = Nt;


%% Input specifications

c = settings.c0;            	%-- speed of sound (m/s)
fs = settings.fs;               %-- sampling frequency (Hz)
theta = settings.PWangle;       %-- PW emission angles (rad)

t = settings.t;                 %-- time: t-axis (s)
z = settings.z;                 %-- imaging depth: z-axis (m)
x = settings.x;                 %-- lateral coordinates: x-axis (m)

dz = mean(diff(z));             %-- depth spacing (m)
dx = mean(diff(x));             %-- lateral spacing (m)

v = c/2;                        %-- ERM velocity assumption (m/s)


%% Configure 1-D FFT grid

NtFFT = settings.NtFFT;         %-- FFT points for t-axis
F0 = (0:NtFFT/2-1)'*fs/NtFFT;	%-- positive temporal frequencies


%% Define slant parameters

pMin = sind(settings.slantMin)/c;	%-- min. slant parameter value
pMax = sind(settings.slantMax)/c;	%-- max. slant parameter value
Np = settings.slantNum;             %-- number of slants

p = linspace(pMin,pMax,Np);         %-- all slant parameter values


%% Initial (zero-valued) beamformed data in z-x domain

BDcomp = zeros(Nz,Nx);              %-- compounded data in z-x domain
BDall = zeros(Nz,Nx,Nf);            %-- all individual frames


%% Slant stacking for each plane wave:

for j = 1:Nf
        
    D_t_x = RFdata(:,:,j);
    
    D_F_x = fft(D_t_x,NtFFT);   
    D_F_x = D_F_x(1:NtFFT/2,:);         %-- extract positive-F spectrum    
    
    D_z_x = zeros(Nz,Nx);    
    for ip = 1:Np,
        %-- note: abs(F0) = F0, as only positive-F spectrum is used
        F_line = F0.*sum(D_F_x.*exp(-2i*pi*F0*x'*p(ip)),2);
        t_line = ifft(F_line,NtFFT,'symmetric');
        t_line = t_line(1:Nt);
        %-- compute z/x-dependent tau components
        ztau = z/c*(1+cos(theta(j)))/2*(1+sqrt(1-p(ip)^2*c^2));
        xtau = x*p(ip);
        %-- perform slant calculations
        slant = zeros(Nz,Nx);
        for ix = 1:Nx,
            tau = ztau + xtau(ix);
            slant(:,ix) = interp1(t,t_line,tau,'linear',0);
        end;
        %-- perform slant stacking
        D_z_x = D_z_x + slant;
    end
    
    %-- compensate for plane-wave angle in z-x domain
	if theta(j)~=0
        for i = 1:Nx
            znew = z + x(i)*tan(theta(j))/2;
            D_z_x(:,i) = interp1(z,D_z_x(:,i),znew,'linear',0);
        end
    end
       
    %-- compound data in z-x domain
    BDcomp = BDcomp + D_z_x;
    BDall(:,:,j) = D_z_x;               %-- optional: save each frame    
    
end

end
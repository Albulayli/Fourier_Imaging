
%--------------------------------------------------------------------------
%-- Demo using PICMUS benchmark data (IEEE IUS 2016): "expres" phantom
%-- Fourier-domain depth migration (FDM): RF data, modified Stolt's method, lense shift compensation 
%-- Authors: M.Albulayli and D.Rakhmatov (daler@ece.uvic.ca)
%-- Institution: ECE Department, University of Victoria
%-- Date: 2018/05/10
%--------------------------------------------------------------------------

clear all;
close all;
clc;

addpath(genpath('./src'));

load('../Reference/expres_scan.mat');	%-- MAT structure identifier: "expres_scan"
rf_scan = expres_scan;
load('../Reference/expres.mat');        %-- MAT structure identifier: "expres"
dataset = expres;


%-- dataset.firings = total number of steered plane-wave emissions
pw_indices{1} = 38;
pw_indices{2} = round(linspace(1,dataset.firings,11));
pw_indices{3} = round(1:dataset.firings);

%-- Configure z-axis (US image depth) and x-axis (ULA sensor position)
t_axis = (0:(size(dataset.data,1)-1)).'/dataset.sampling_frequency+dataset.initial_time;
z_axis = t_axis*dataset.c0/2;
x_axis = dataset.probe_geometry(:,1);
rx_f_number = 1.75;


%% Define scan structure

scan_x_axis = x_axis;
scan_z_axis = z_axis((min(rf_scan.z_axis) <= z_axis) & (z_axis <= max(rf_scan.z_axis)));
[scan_x_matrix,scan_z_matrix] = meshgrid(scan_x_axis,scan_z_axis);
scan_x = scan_x_matrix(:);
scan_z = scan_z_matrix(:);
scan_pixels = length(scan_x);

scan.x_axis = scan_x_axis;
scan.z_axis = scan_z_axis;
scan.x_matrix = scan_x_matrix;
scan.z_matrix = scan_z_matrix;
scan.x = scan_x;
scan.z = scan_z;
scan.dx = mean(diff(scan_x_axis));
scan.dz = mean(diff(scan_z_axis));
scan.Nx = numel(scan_x_axis);
scan.Nz = numel(scan_z_axis);
scan.pixels = scan_pixels;


%% Configure settings and perform beamforming

settings.fs = dataset.sampling_frequency;
settings.c0 = dataset.c0;
settings.t = t_axis;
settings.z = z_axis;
settings.x = x_axis;
settings.NtFFT = pow2(nextpow2(length(t_axis)+1));      %-- number of temporal FFT points
settings.NxFFT = pow2(nextpow2(length(x_axis)+1));      %-- number of spatial FFT points

beamformed_data = zeros(length(z_axis),length(x_axis),length(pw_indices));
envelope_beamformed_data = zeros(length(z_axis),length(x_axis),length(pw_indices));

for f=1:length(pw_indices)
    pw = pw_indices{f};
	RFdata = dataset.data(:,:,pw);
    settings.PWangle = dataset.angles(pw);
    lense_shift = 32;   %-- lense_shift/probe.fs*probe.c0 = 2.3656 mm (PICMUS)
    for k=1:length(pw)
        RFdata(:,:,k) = [zeros(lense_shift,size(RFdata,2)); RFdata(1:end-lense_shift,:,k)];
    end 
    [BDcomp,~] = pwStolt(RFdata,settings);
    beamformed_data(:,:,f) = [BDcomp(lense_shift+1:end,:); zeros(lense_shift,size(beamformed_data,2))];
    envelope_beamformed_data(:,:,f) = abs(hilbert(beamformed_data(:,:,f)));
    disp([num2str(pw),' / ',num2str(length(pw))])
end

save('expres_bfdata_pwStolt.mat', 'beamformed_data');


%% Interpolate envelope from [z_axis,x_axis] into [scan.z_axis,scan.x_axis]

%-- interpolate the requested grid
resampled_envelope_beamformed_data = zeros(numel(scan.z_axis),numel(scan.x_axis),numel(pw_indices));

[Z,X] = ndgrid(z_axis,x_axis);
[iZ,iX] = ndgrid(scan.z_axis,scan.x_axis);
for f=1:length(pw_indices)
    gI = griddedInterpolant(Z,X,envelope_beamformed_data(:,:,f),'linear');
    resampled_envelope_beamformed_data(:,:,f) = gI(iZ,iX);
end



%% Save image structure for evaluation

expres_image.scan = scan;
expres_image.number_plane_waves = cellfun('length',pw_indices);
expres_image.data = resampled_envelope_beamformed_data;
expres_image.transmit_f_number = 0;
expres_image.receive_f_number = rx_f_number;
expres_image.method = 'PW Stolt';

save('expres_image_pwStolt.mat', 'expres_image');


%% Display interpolated and original images (60-dB range)

vrange = [-60 0];
x_lim = [min(scan.x_axis) max(scan.x_axis)]*1e3;
z_lim = [min(scan.z_axis) max(scan.z_axis)]*1e3;

% %-- show original images
% for f=1:length(pw_indices)    
%     %-- compute dB values
%     env = envelope_beamformed_data(:,:,f);
%     im = 20*log10(env./max(env(:)));     
%     %-- display image    
%     figure(f);
%     imagesc((x_axis)*1e3,(z_axis)*1e3,im);
%     shading flat; colormap gray; caxis(vrange); %colorbar;
%     set(gca,'FontName','Courier','FontSize',8,'FontWeight','demi');
%     axis tight image;
%     xlabel('x (mm)','FontName','Courier','FontSize',10,'FontWeight','demi');
%     ylabel('z (mm)','FontName','Courier','FontSize',10,'FontWeight','demi');
%     set(gca,'YDir','reverse');
%     npw = length(pw_indices{f});
%     if npw == 1, title('Stolt (Ours), 1 Plane Wave','FontName','Courier','FontSize',10,'FontWeight','bold');
%     else title(sprintf('Stolt (Ours), %d Plane Waves',npw),'FontName','Courier','FontSize',10,'FontWeight','bold');
%     end;
%     saveas(gcf,sprintf('expres_pwStolt_original_%d.fig',npw));
%     pause(0.5);
% end

%-- show interpolated images
for f=1:length(pw_indices)    
    %-- compute dB values
    env = resampled_envelope_beamformed_data(:,:,f);
    im = 20*log10(env./max(env(:)));   
    %-- display image    
    figure(f);
    imagesc((scan.x_axis)*1e3,(scan.z_axis)*1e3,im);
    shading flat; colormap gray; caxis(vrange); colorbar;
    set(gca,'FontName','Courier','FontSize',16,'FontWeight','demi');
    axis([x_lim z_lim]);
    axis equal manual;
    xlabel('x (mm)','FontName','Courier','FontSize',16,'FontWeight','demi');
    ylabel('z (mm)','FontName','Courier','FontSize',16,'FontWeight','demi');
    set(gca,'YDir','reverse');
    npw = length(pw_indices{f});
    if npw == 1, title('Stolt (Ours), One PW','FontName','Courier','FontSize',16,'FontWeight','bold');
    else title(sprintf('Stolt (Ours), %d PWs',npw),'FontName','Courier','FontSize',16,'FontWeight','bold');
    end;
    saveas(gcf,sprintf('expres_pwStolt_%d.fig',npw));
    saveas(gcf,sprintf('expres_pwStolt_%d',npw),'eps');
    pause(0.5);
end


%% Perform image evaluation

addpath(genpath('../Reference'));
exec_evaluation_resolution_distorsion_exp('expres_image_pwStolt.mat','expres_result_pwStolt.txt');


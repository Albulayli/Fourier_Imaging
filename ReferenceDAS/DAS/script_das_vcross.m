%-- Function which implements the conventional Delay And Sum (DAS) beamform technique with apodization in reception
%-- The corresponding code is dedicated to the reconstrucion of dataset (rawdata) saved in RF format

%-- Authors: Alfonso Rodriguez-Molares (alfonso.r.molares@ntnu.no)
%--          Olivier Bernard (olivier.bernard@creatis.insa-lyon.fr)

%-- $Date: 2016/03/01 $


%-- Modified on 2018/05/11 by D.N.Rakhmatov (daler@ece.uvic.ca)
%-- DAS beamforming: RF data, no windowing, lense shift compensation, "vcross" phantom

clear all;
close all;
clc;

addpath(genpath('../Reference'));

load('../Reference/vcross_scan.mat');	%-- MAT structure identifier: "vcross_scan"
rf_scan = vcross_scan;
load('../Reference/vcross.mat');        %-- MAT structure identifier: "vcross"
dataset = vcross;


%     assert(isempty(dataset.modulation_frequency)||dataset.modulation_frequency==0,'The supplied dataset is not RF');
%
%     %-- select the plane waves that will be used in each frame
%     if nargin < 3
%         pw_indices{1} = 1:dataset.firings;
%     end

pw_indices{1} = 38;
pw_indices{2} = round(linspace(1,dataset.firings,11));
pw_indices{3} = round(1:dataset.firings);               %-- dataset.firings corresponding to the total number of emitted steered plane waves

%-- define scan based on time axis
time = (0:(size(dataset.data,1)-1)).'/dataset.sampling_frequency+dataset.initial_time;
z_axis = time*dataset.c0/2;
%     rf_scan = linear_scan(scan.x_axis,z_axis);

x_axis = dataset.probe_geometry(:,1);
[x_matrix,z_matrix] = meshgrid(x_axis,z_axis);
x = x_matrix(:);
z = z_matrix(:);
pixels = length(x);


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


%% Perform DAS beamforming

%-- receive apodization
%-- dynamically expanding receive aperture with hanning apodization
rx_f_number = 1.75;
%     rx_aperture = rf_scan.z/rx_f_number;
rx_aperture = z/rx_f_number;
%     rx_aperture_distance = abs(rf_scan.x*ones(1,dataset.channels)-ones(rf_scan.pixels,1)*dataset.probe_geometry(:,1).');
rx_aperture_distance = abs(x*ones(1,dataset.channels)-ones(pixels,1)*dataset.probe_geometry(:,1).');

%     receive_apodization = tools.apodization(rx_aperture_distance,rx_aperture*ones(1,dataset.channels),'hanning');
receive_apodization = tools.apodization(rx_aperture_distance,rx_aperture*ones(1,dataset.channels),'hamming');

%-- angular apodization -> no apodization
%     angular_apodization = ones(rf_scan.pixels,dataset.firings);
angular_apodization = ones(pixels,dataset.firings);

%-- beamforming loop
%     beamformed_data = zeros(rf_scan.pixels,length(pw_indices));
beamformed_data = zeros(pixels,length(pw_indices));
time_vector = dataset.initial_time+(0:(dataset.samples-1))/dataset.sampling_frequency;
%     wb=waitbar(0,'DAS beamforming');


lense_shift = 32;	%-- lense_shift/probe.fs*probe.c0 = 2.3656 mm (PICMUS)

for f=1:length(pw_indices)
    %         waitbar(f/length(pw_indices),wb,sprintf('DAS-RF beamforming %0.0f%%',f/length(pw_indices)*100));
    for pw=pw_indices{f}
        
        rf_data = dataset.data(:,:,pw);
        %-- compensate for lense shift
        rf_data = [zeros(lense_shift,size(rf_data,2)); rf_data(1:end-lense_shift,:)];
        
        %-- transmit delay
%         transmit_delay = rf_scan.z*cos(dataset.angles(pw))+rf_scan.x*sin(dataset.angles(pw));
        transmit_delay = z*cos(dataset.angles(pw))+x*sin(dataset.angles(pw));        
        for nrx=1:dataset.channels
            %-- receive delay
%             receive_delay = sqrt((dataset.probe_geometry(nrx,1)-rf_scan.x).^2+(dataset.probe_geometry(nrx,3)-rf_scan.z).^2);
            receive_delay = sqrt((dataset.probe_geometry(nrx,1)-x).^2+(dataset.probe_geometry(nrx,3)-z).^2);            
            %-- total delay
            delay = (transmit_delay+receive_delay)/dataset.c0;
            %-- beamformed data
%             beamformed_data(:,f) = beamformed_data(:,f)+angular_apodization(:,pw).*receive_apodization(:,nrx).*interp1(time_vector,dataset.data(:,nrx,pw),delay,'spline',0);
            beamformed_data(:,f) = beamformed_data(:,f)+angular_apodization(:,pw).*receive_apodization(:,nrx).*interp1(time_vector,rf_data(:,nrx),delay,'spline',0);
        end
        %             clc;
        disp([num2str(pw),' / ',num2str(length(pw_indices{f}))])
    end
end
%     close(wb);

beamformed_data(isnan(beamformed_data))=0;

%-- reshape
%     reshaped_beamformed_data = reshape(beamformed_data,[numel(rf_scan.z_axis) numel(rf_scan.x_axis)  length(pw_indices)]);
reshaped_beamformed_data = reshape(beamformed_data,[numel(z_axis) numel(x_axis) length(pw_indices)]);


for f=1:length(pw_indices)
    %-- compensate for lense shift
    reshaped_beamformed_data(:,:,f) = [reshaped_beamformed_data(lense_shift+1:end,:,f); zeros(lense_shift,size(reshaped_beamformed_data,2))];
end


save('vcross_bfdata_das.mat', 'reshaped_beamformed_data');


%-- compute envelope
%     envelope_beamformed_data = tools.envelope(reshaped_beamformed_data);
envelope_beamformed_data = reshape(abs(hilbert(reshaped_beamformed_data(:,:))),size(reshaped_beamformed_data));


%% Interpolate envelope from [z_axis,x_axis] into [scan.z_axis,scan.x_axis]

%-- interpolate the requested grid
resampled_envelope_beamformed_data = zeros(numel(scan.z_axis),numel(scan.x_axis),numel(pw_indices));

% for f=1:length(pw_indices)
%     resampled_envelope_beamformed_data(:,:,f) = interp1(rf_scan.z_axis,envelope_beamformed_data(:,:,f),scan.z_axis,'linear',0);
% end

[Z,X] = ndgrid(z_axis,x_axis);
[iZ,iX] = ndgrid(scan.z_axis,scan.x_axis);
for f=1:length(pw_indices)
    gI = griddedInterpolant(Z,X,envelope_beamformed_data(:,:,f),'linear');
    resampled_envelope_beamformed_data(:,:,f) = gI(iZ,iX);
end


%     %-- declare an us_image object to store the beamformed data
%     image = us_image('DAS-RF beamforming');
%     image.author = 'Alfonso Rodriguez-Molares <alfonso.r.molares@ntnu.no>';
%     image.affiliation = 'Norwegian University of Science and Technology (NTNU)';
%     image.algorithm = 'Delay-and-Sum (RF version)';
%     image.scan = scan;
%     image.number_plane_waves = cellfun('length',pw_indices);
%     image.data = resampled_envelope_beamformed_data;
%     image.transmit_f_number = 0;
%     image.receive_f_number = rx_f_number;
%     image.transmit_apodization_window = 'none';
% %     image.receive_apodization_window = 'Tukey 25%';
%     image.receive_apodization_window = 'none';



%% Save image structure for evaluation

vcross_image.scan = scan;
vcross_image.number_plane_waves = cellfun('length',pw_indices);
vcross_image.data = resampled_envelope_beamformed_data;
vcross_image.transmit_f_number = 0;
vcross_image.receive_f_number = rx_f_number;
vcross_image.method = 'DAS beamforming (RF data, no receive apodization)';

save('vcross_image_das.mat', 'vcross_image');


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
%     if npw == 1, title('DAS, 1 Plane Wave','FontName','Courier','FontSize',10,'FontWeight','bold');
%     else title(sprintf('DAS, %d Plane Waves',npw),'FontName','Courier','FontSize',10,'FontWeight','bold');
%     end;
%     saveas(gcf,sprintf('vcross_das_original_%d.fig',npw));
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
    if npw == 1, title('DAS, One PW','FontName','Courier','FontSize',16,'FontWeight','bold');
    else title(sprintf('DAS, %d PWs',npw),'FontName','Courier','FontSize',16,'FontWeight','bold');
    end;
    saveas(gcf,sprintf('vcross_das_%d.fig',npw));
    saveas(gcf,sprintf('vcross_das_%d',npw),'eps');
    pause(0.5);
end

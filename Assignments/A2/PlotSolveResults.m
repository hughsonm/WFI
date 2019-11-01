clearvars;
close all;
%% Read in Data
ITX = 1;
FREQ = 1E9;
CNAUGHT = 2.9979E8;
K_B = (2*pi*FREQ)/CNAUGHT;

LAMBDA = CNAUGHT/FREQ;

%% Plot Forward Problem
FWD_DIR = 'FwdData_SingleTx/';
fwd_eps_r =     ReadCppVectorFromFile([FWD_DIR 'eps_r.txt']);
fwd_k2 =        fwd_eps_r*(K_B)^2;
fwd_tri =       ReadCppMatrixFromFile([FWD_DIR 'tri_tri.txt']);
fwd_pts =       ReadCppMatrixFromFile([FWD_DIR 'tri_pts.txt']);

fwd_tags = fwd_tri(:,end);
fwd_tri = fwd_tri(:,1:end-1);
fwd_tri = fwd_tri+1;

% eq_radii = sqrt(areas/pi);

Ez_inc =    ReadCppMatrixFromFile([FWD_DIR 'Ez_inc.txt']);
Ez_sct =    ReadCppMatrixFromFile([FWD_DIR 'Ez_sct.txt']);
Ez_tot =    ReadCppMatrixFromFile([FWD_DIR 'Ez_tot.txt']);



% probes =    ReadCppMatrixFromFile('ProbePositions4A2.txt');
Ez_inc_p =  ReadCppMatrixFromFile([FWD_DIR 'Ez_inc_d.txt']);
Ez_sct_p =  ReadCppMatrixFromFile([FWD_DIR 'Ez_sct_d.txt']);
Ez_tot_p =  ReadCppMatrixFromFile([FWD_DIR 'Ez_tot_d.txt']);

xyz_p =     ReadCppMatrixFromFile([FWD_DIR 'probe_xyz.txt']);

fwd_fig = figure();

subplot(4,2,1);
trisurf(fwd_tri,fwd_pts(:,1),fwd_pts(:,2),fwd_pts(:,3),real(fwd_eps_r));
title('Re[\epsilon_r]');
xlabel('x [m]');ylabel('y [m]');zlabel('z [m]');
colorbar;axis image; view(2);

subplot(4,2,2);
trisurf(fwd_tri,fwd_pts(:,1),fwd_pts(:,2),fwd_pts(:,3),imag(fwd_eps_r));
title('Im[\epsilon_r]');
xlabel('x [m]');ylabel('y [m]');zlabel('z [m]');
colorbar;axis image; view(2);

subplot(4,2,3);
trisurf(fwd_tri,fwd_pts(:,1),fwd_pts(:,2),fwd_pts(:,3),real(Ez_inc(:,ITX)),...
    'LineStyle','None');
title('Re[E_z^{inc}]');
xlabel('x [m]');ylabel('y [m]');zlabel('z [m]');
colorbar;
axis image;
view(2);

subplot(4,2,4);
trisurf(fwd_tri,fwd_pts(:,1),fwd_pts(:,2),fwd_pts(:,3),imag(Ez_inc(:,ITX)),...
    'LineStyle','None');
title('Im[E_z^{inc}]');
xlabel('x [m]');ylabel('y [m]');zlabel('z [m]');
colorbar;
axis image;
view(2);

subplot(4,2,5);
trisurf(fwd_tri,fwd_pts(:,1),fwd_pts(:,2),fwd_pts(:,3),real(Ez_tot(:,ITX)),...
    'LineStyle','None');
title('Re[E_z^{tot}]');
xlabel('x [m]');ylabel('y [m]');zlabel('z [m]');
colorbar;
axis image;
view(2);

subplot(4,2,6);
trisurf(fwd_tri,fwd_pts(:,1),fwd_pts(:,2),fwd_pts(:,3),imag(Ez_tot(:,ITX)),...
    'LineStyle','None');
title('Im[E_z^{tot}]');
xlabel('x [m]');ylabel('y [m]');zlabel('z [m]');
colorbar;
axis image;
view(2);

subplot(4,2,7);
trisurf(fwd_tri,fwd_pts(:,1),fwd_pts(:,2),fwd_pts(:,3),real(Ez_sct(:,ITX)),...
    'LineStyle','None');
title('Re[E_z^{sct}]');
xlabel('x [m]');ylabel('y [m]');zlabel('z [m]');
colorbar;
axis image;
view(2);

subplot(4,2,8);
trisurf(fwd_tri,fwd_pts(:,1),fwd_pts(:,2),fwd_pts(:,3),imag(Ez_sct(:,ITX)),...
    'LineStyle','None');
title('Im[E_z^{sct}]');
xlabel('x [m]');ylabel('y [m]');zlabel('z [m]');
colorbar;
axis image;
view(2);


probe_fig = figure();
sep = 2*max(max(abs(Ez_inc_p)));
subplot(3,2,1);
for ii = 1:size(Ez_inc_p,2)
    plot3(xyz_p(:,1),xyz_p(:,2),real(Ez_inc_p(:,ii))+ii*sep);hold on;
end
title('Re[E_z^{inc}] at probes');
grid on;

subplot(3,2,2);
for ii = 1:size(Ez_inc_p,2)
    plot3(xyz_p(:,1),xyz_p(:,2),imag(Ez_inc_p(:,ii))+ii*sep);hold on;
end
title('Im[E_z^{inc}] at probes');
grid on;

sep = 2*max(max(abs(Ez_tot_p)));
subplot(3,2,3);
for ii = 1:size(Ez_tot_p,2)
    plot3(xyz_p(:,1),xyz_p(:,2),real(Ez_tot_p(:,ii))+ii*sep);hold on;
end
title('Re[E_z^{tot}] at probes');
grid on;

subplot(3,2,4);
for ii = 1:size(Ez_tot_p,2)
    plot3(xyz_p(:,1),xyz_p(:,2),imag(Ez_tot_p(:,ii))+ii*sep);hold on;
end
title('Im[E_z^{tot}] at probes');
grid on;

sep = 2*max(max(abs(Ez_sct_p)));
subplot(3,2,5);
for ii = 1:size(Ez_sct_p,2)
    plot3(xyz_p(:,1),xyz_p(:,2),real(Ez_sct_p(:,ii))+ii*sep);hold on;
end
title('Re[E_z^{sct}] at probes');
grid on;

subplot(3,2,6);
for ii = 1:size(Ez_sct_p,2)
    plot3(xyz_p(:,1),xyz_p(:,2),imag(Ez_sct_p(:,ii))+ii*sep);hold on;
end
title('Im[E_z^{sct}] at probes');
grid on;


%% Plot Inversion
INV_DIR = 'A2Output_MultiTx/';
inv_tri =       ReadCppMatrixFromFile([INV_DIR 'Tri.txt']);
inv_pts =       ReadCppMatrixFromFile([INV_DIR 'Points.txt']);
q5_w_tikh =     ReadCppMatrixFromFile([INV_DIR 'q5_w_tikh.txt']);
q5_X_tikh =     ReadCppVectorFromFile([INV_DIR 'q5_X_tikh.txt']);
Ez_opt_tikh =   ReadCppMatrixFromFile([INV_DIR 'q5_u_tikh.txt']);
q5_w_unreg =    ReadCppMatrixFromFile([INV_DIR 'q5_w_unreg.txt']);
q5_X_unreg =    ReadCppVectorFromFile([INV_DIR 'q5_X_unreg.txt']);
Ez_opt_unreg =  ReadCppMatrixFromFile([INV_DIR 'q5_u_unreg.txt']);

q3_w =          ReadCppMatrixFromFile([INV_DIR 'q3_w.txt']);
q3_X =          ReadCppVectorFromFile([INV_DIR 'q3_X.txt']);

Gbd =           ReadCppMatrixFromFile([INV_DIR 'Gbd.txt']);


inv_tags = inv_tri(:,end);
inv_tri = inv_tri(:,1:end-1);
inv_tri = inv_tri+1; % Add one because this was made by a zero-based indexing code

ann_w_fig = figure();
subplot(2,1,1);
trisurf(inv_tri,...
    inv_pts(:,1),inv_pts(:,2),0*inv_pts(:,1),...
    real(q3_w(:,1)));
title('Re[w] at tx 1 by annihilation');
axis image; view(2);colorbar;
subplot(2,1,2);
trisurf(inv_tri,...
    inv_pts(:,1),inv_pts(:,2),0*inv_pts(:,1),...
    imag(q3_w(:,1)));
title('Im[w] at tx 1 by annihilation');
axis image; view(2);colorbar;

ann_X_fig = figure();
subplot(2,1,1);
trisurf(inv_tri,...
    inv_pts(:,1),inv_pts(:,2),0*inv_pts(:,1),...
    db(abs(real(q3_X))));
title('Real part of contrast by annihilation');
axis image; view(2);colorbar;
subplot(2,1,2);
trisurf(inv_tri,...
    inv_pts(:,1),inv_pts(:,2),0*inv_pts(:,1),...
    db(abs(imag(q3_X))));
title('Imag part of contrast by annihilation');
axis image; view(2);colorbar;

pb_unreg_figs = [];
for itx = 1:4:size(Ez_opt_tikh,2)
%     itx = 1;
    pb_unreg_figs(end+1) = figure();
    subplot(4,1,1);
    
    trisurf(inv_tri,inv_pts(:,1),inv_pts(:,2),0*inv_pts(:,1),real(q5_w_unreg(:,itx))); view(2);colorbar;axis image;
    title(['Re[w] for tx ' num2str(itx)]);     
    
    subplot(4,1,2);
    
    trisurf(inv_tri,inv_pts(:,1),inv_pts(:,2),0*inv_pts(:,1),imag(q5_w_unreg(:,itx))); view(2);colorbar;axis image;
    title(['Im[w] for tx ' num2str(itx)]);            
    
    subplot(4,1,3);
    trisurf(inv_tri,inv_pts(:,1),inv_pts(:,2),0*inv_pts(:,1),real(q5_X_unreg));view(2);colorbar;axis image;
    title('Re[\chi]');
    
    subplot(4,1,4);
    trisurf(inv_tri,inv_pts(:,1),inv_pts(:,2),0*inv_pts(:,1),imag(q5_X_unreg));view(2);colorbar;axis image;
    title('Im[\chi]');
end

pb_tikh_figs = [];
for itx = 1:4:size(Ez_opt_tikh,2)
%     itx = 1;
    pb_tikh_figs(end+1) = figure();
    subplot(4,1,1);
    
    trisurf(inv_tri,inv_pts(:,1),inv_pts(:,2),0*inv_pts(:,1),real(q5_w_tikh(:,itx))); view(2);colorbar;axis image;
    title(['Re[w] for tx ' num2str(itx)]);     
    
    subplot(4,1,2);
    
    trisurf(inv_tri,inv_pts(:,1),inv_pts(:,2),0*inv_pts(:,1),imag(q5_w_tikh(:,itx))); view(2);colorbar;axis image;
    title(['Im[w] for tx ' num2str(itx)]);            
    
    subplot(4,1,3);
    trisurf(inv_tri,inv_pts(:,1),inv_pts(:,2),0*inv_pts(:,1),real(q5_X_tikh));view(2);colorbar;axis image;
    title('Re[\chi]');
    
    subplot(4,1,4);
    trisurf(inv_tri,inv_pts(:,1),inv_pts(:,2),0*inv_pts(:,1),imag(q5_X_tikh));view(2);colorbar;axis image;
    title('Im[\chi]');
end

curve_fig = figure();
ncurves = size(q5_w_tikh,2);
for itx = 1:ncurves
    cc = ReadCppMatrixFromFile([INV_DIR 'q5_curve_' num2str(itx-1) '.txt']);
    pp = plot(cc(1,:),cc(2,:),'LineWidth',5);
    pp.Color(4) = 1/ncurves;
    np = size(cc,2);
    minexp = -floor(np/2);
    maxexp = floor(np/2);
    text(cc(1,:),cc(2,:),num2cell(10.^(minexp:maxexp)));
    hold on;
    xlabel('log(|d-P*w|)');
    ylabel('log(|w|)');
end
title('Tikhonov Curves For All Transmitters');

if(1 < size(q5_w_tikh,2))
    fig_dir = './Report/figs/Multi/';
else
    fig_dir = './Report/figs/Single/';
end
mkdir(fig_dir);
fig_ext = '.png';
drawnow;
figure(fwd_fig);saveas(fwd_fig,[fig_dir 'Forward' fig_ext]);
figure(probe_fig);saveas(probe_fig,[fig_dir 'Probe' fig_ext]);
figure(ann_w_fig);saveas(ann_w_fig,[fig_dir 'AnnW' fig_ext]);
figure(ann_X_fig);saveas(ann_X_fig,[fig_dir 'AnnX' fig_ext]);
for ii = 1:length(pb_unreg_figs)
    figure(pb_unreg_figs(ii));saveas(pb_unreg_figs(ii),[fig_dir 'PB_unreg' num2str(ii) fig_ext]);
    figure(pb_tikh_figs(ii));saveas(pb_tikh_figs(ii),[fig_dir 'PB_tikh' num2str(ii) fig_ext]);
end
figure(curve_fig);saveas(curve_fig,[fig_dir 'Curve' fig_ext]);






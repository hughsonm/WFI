clearvars;
close all;
%% Read in Data
ITX = 1;
FREQ = 5E9;
CNAUGHT = 2.9979E8;
K_B = (2*pi*FREQ)/CNAUGHT;

LAMBDA = CNAUGHT/FREQ;

%% Plot Forward Problem
fwd_eps_r =     ReadCppVectorFromFile('eps_r.txt');
fwd_k2 =        fwd_eps_r*(K_B)^2;
fwd_tri =       ReadCppMatrixFromFile('tri_tri.txt');
fwd_pts =       ReadCppMatrixFromFile('tri_pts.txt');

fwd_tags = fwd_tri(:,end);
fwd_tri = fwd_tri(:,1:end-1);
fwd_tri = fwd_tri+1;

% eq_radii = sqrt(areas/pi);

Ez_inc =    ReadCppMatrixFromFile('Ez_inc.txt');
Ez_sct =    ReadCppMatrixFromFile('Ez_sct.txt');
Ez_tot =    ReadCppMatrixFromFile('Ez_tot.txt');



% probes =    ReadCppMatrixFromFile('ProbePositions4A2.txt');
Ez_inc_p =  ReadCppMatrixFromFile('Ez_inc_d.txt');
Ez_sct_p =  ReadCppMatrixFromFile('Ez_sct_d.txt');
Ez_tot_p =  ReadCppMatrixFromFile('Ez_tot_d.txt');

xyz_p =     ReadCppMatrixFromFile('probe_xyz.txt');

figure();

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


figure();
sep = 2*max(max(abs(Ez_inc_p)));
subplot(3,2,1);
for ii = 1:size(Ez_inc_p,2)
    plot3(xyz_p(:,1),xyz_p(:,2),real(Ez_inc_p(:,ii))+ii*sep);hold on;
end
title('Re[E_z^{inc}] at probes');
subplot(3,2,2);
for ii = 1:size(Ez_inc_p,2)
    plot3(xyz_p(:,1),xyz_p(:,2),imag(Ez_inc_p(:,ii))+ii*sep);hold on;
end
title('Im[E_z^{inc}] at probes');

sep = 2*max(max(abs(Ez_tot_p)));
subplot(3,2,3);
for ii = 1:size(Ez_tot_p,2)
    plot3(xyz_p(:,1),xyz_p(:,2),real(Ez_tot_p(:,ii))+ii*sep);hold on;
end
title('Re[E_z^{tot}] at probes');
subplot(3,2,4);
for ii = 1:size(Ez_tot_p,2)
    plot3(xyz_p(:,1),xyz_p(:,2),imag(Ez_tot_p(:,ii))+ii*sep);hold on;
end
title('Im[E_z^{tot}] at probes');

sep = 2*max(max(abs(Ez_sct_p)));
subplot(3,2,5);
for ii = 1:size(Ez_sct_p,2)
    plot3(xyz_p(:,1),xyz_p(:,2),real(Ez_sct_p(:,ii))+ii*sep);hold on;
end
title('Re[E_z^{sct}] at probes');
subplot(3,2,6);
for ii = 1:size(Ez_sct_p,2)
    plot3(xyz_p(:,1),xyz_p(:,2),imag(Ez_sct_p(:,ii))+ii*sep);hold on;
end
title('Im[E_z^{sct}] at probes');




%% Plot Inversion
inv_tri =       ReadCppMatrixFromFile('Tri.txt');
inv_pts =       ReadCppMatrixFromFile('Points.txt');
alphas =    ReadCppMatrixFromFile('Alphas.txt');
X =         ReadCppVectorFromFile('Chi.txt');
Ez_opt =    ReadCppMatrixFromFile('Ez_tot_opt.txt');
inv_centroids = ReadCppMatrixFromFile('tri_centroids.txt');
inv_areas =     ReadCppVectorFromFile('tri_areas.txt');


inv_tags = inv_tri(:,end);
inv_tri = inv_tri(:,1:end-1);
inv_tri = inv_tri+1; % Add one because this was made by a zero-based indexing code



for itx = 1:4:size(Ez_opt,2)
%     itx = 1;
    figure();
    subplot(4,2,1);
    
    trisurf(inv_tri,inv_pts(:,1),inv_pts(:,2),0*inv_pts(:,1),real(alphas(:,itx))); view(2);colorbar;axis image;
    title(['Re[w] for tx ' num2str(itx)]);     
    
    subplot(4,2,2);
    
    trisurf(inv_tri,inv_pts(:,1),inv_pts(:,2),0*inv_pts(:,1),imag(alphas(:,itx))); view(2);colorbar;axis image;
    title(['Im[w] for tx ' num2str(itx)]);            
    
    subplot(4,2,3);
    trisurf(inv_tri,inv_pts(:,1),inv_pts(:,2),0*inv_pts(:,1),real(X));view(2);colorbar;axis image;
    title('Re[\chi]');
    
    subplot(4,2,4);
    trisurf(inv_tri,inv_pts(:,1),inv_pts(:,2),0*inv_pts(:,1),imag(X));view(2);colorbar;axis image;
    title('Im[\chi]');
            
    subplot(4,2,5);
    trisurf(inv_tri,inv_pts(:,1),inv_pts(:,2),0*inv_pts(:,1),real(Ez_opt(:,itx)));view(2);colorbar;axis image;
    title(['Optimal Re[E_{z}^{tot}] for tx ' num2str(itx)]);
    
    subplot(4,2,6);
    trisurf(inv_tri,inv_pts(:,1),inv_pts(:,2),0*inv_pts(:,1),imag(Ez_opt(:,itx)));view(2);colorbar;axis image;
    title(['Optimal Im[E_{z}^{tot}] for tx ' num2str(itx)]);
    
    subplot(4,2,7);
    trisurf(inv_tri,inv_pts(:,1),inv_pts(:,2),0*inv_pts(:,1),real(alphas(:,itx)./Ez_opt(:,itx)));view(2);colorbar;axis image;
    title(['Re[W/E] for tx ' num2str(itx)]);
    
    subplot(4,2,8);
    trisurf(inv_tri,inv_pts(:,1),inv_pts(:,2),0*inv_pts(:,1),imag(alphas(:,itx)./Ez_opt(:,itx)));view(2);colorbar;axis image;
    title(['Im[W/E] for tx ' num2str(itx)]);
end






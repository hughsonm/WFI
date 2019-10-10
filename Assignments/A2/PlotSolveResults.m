clearvars;
close all;
%% Read in Data
ITX = 1;
FREQ = 5E9;
CNAUGHT = 2.9979E8;
K_B = (2*pi*FREQ)/CNAUGHT;

LAMBDA = CNAUGHT/FREQ;

%% Plot Forward Problem
fwd_k2 =        ReadCppVectorFromFile('k2_fgd.txt');
fwd_eps_r =     fwd_k2/(K_B)^2;
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

figure();

subplot(2,2,1);
trisurf(fwd_tri,fwd_pts(:,1),fwd_pts(:,2),fwd_pts(:,3),real(fwd_eps_r));
title('eps_r');
xlabel('x [m]');ylabel('y [m]');zlabel('z [m]');
colorbar;axis image; view(2);


subplot(2,2,2);
trisurf(fwd_tri,fwd_pts(:,1),fwd_pts(:,2),fwd_pts(:,3),abs(Ez_inc(:,ITX)),...
    'LineStyle','None');
title('|E_z^{inc}|');
xlabel('x [m]');ylabel('y [m]');zlabel('z [m]');
colorbar;
axis image;
view(2);

subplot(2,2,3);
trisurf(fwd_tri,fwd_pts(:,1),fwd_pts(:,2),fwd_pts(:,3),abs(Ez_sct(:,ITX)),...
    'LineStyle','None');
title('|E_z^{sct}|');
xlabel('x [m]');ylabel('y [m]');zlabel('z [m]');
colorbar;
axis image;
view(2);

subplot(2,2,4);
trisurf(fwd_tri,fwd_pts(:,1),fwd_pts(:,2),fwd_pts(:,3),abs(Ez_tot(:,ITX)),...
    'LineStyle','None');
title('|E_z^{tot}|');
xlabel('x [m]');ylabel('y [m]');zlabel('z [m]');
colorbar;
axis image;
view(2);


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



for itx = 1:12:size(Ez_opt,2)
%     itx = 1;
    figure();
    subplot(3,2,1);
    
    trisurf(inv_tri,inv_pts(:,1),inv_pts(:,2),0*inv_pts(:,1),real(alphas(:,itx))); view(2);colorbar;axis image;
    title(['Contrast Sources for tx ' num2str(itx)]);
    
    subplot(3,2,3);
    trisurf(inv_tri,inv_pts(:,1),inv_pts(:,2),0*inv_pts(:,1),real(X));view(2);colorbar;axis image;
    title(['Real \chi for tx ' num2str(itx)]);
    
    subplot(3,2,4);
    trisurf(inv_tri,inv_pts(:,1),inv_pts(:,2),0*inv_pts(:,1),imag(X));view(2);colorbar;axis image;
    title(['Imag \chi for tx ' num2str(itx)]);
    
    subplot(3,2,5);
    trisurf(inv_tri,inv_pts(:,1),inv_pts(:,2),0*inv_pts(:,1),real(Ez_opt(:,itx)));view(2);colorbar;axis image;
    title(['Optimal E_{z}^{tot} for tx ' num2str(itx)]);
    
    subplot(3,2,6);
    trisurf(inv_tri,inv_pts(:,1),inv_pts(:,2),0*inv_pts(:,1),real(alphas(:,itx)./Ez_opt(:,itx)));view(2);colorbar;axis image;
    title(['W/E for tx ' num2str(itx)]);
end






clearvars;

ITX = 1;
FREQ = 1E9;
CNAUGHT = 2.9979E8;
K_B = (2*pi*FREQ)/CNAUGHT;

LAMBDA = CNAUGHT/FREQ;

tri = dlmread('Test_tri.txt');
pts = dlmread('Test_pts.txt');
centroids = dlmread('tri_centroids.txt');
areas = dlmread('tri_areas.txt');
eq_radii = sqrt(areas/pi);

Ez_inc = dlmread('Ez_inc.txt');
Ez_sct = dlmread('Ez_sct.txt');
Ez_tot = dlmread('Ez_tot.txt');

k2 = dlmread('k2_fgd.txt');
eps_r = k2/K_B;

probes = dlmread('probe_xyz.txt');
Ez_inc_p = dlmread('Ez_inc_d.txt');
Ez_sct_p = dlmread('Ez_sct_d.txt');
Ez_tot_p = dlmread('Ez_tot_d.txt');


tags = tri(:,end);
tri = tri(:,1:end-1);
tri = tri+1; % Add one because this was made by a zero-based indexing code

figure(1)
subplot(2,2,1);
trisurf(tri,pts(:,1),pts(:,2),pts(:,3),abs(k2));
title('k^2');
xlabel('x [m]');ylabel('y [m]');zlabel('z [m]');
colorbar;axis image; view(2);

subplot(2,2,2);
trisurf(tri,pts(:,1),pts(:,2),pts(:,3),abs(Ez_inc(:,ITX)),...
    'LineStyle','None');
title('E_z^{inc}');
xlabel('x [m]');ylabel('y [m]');zlabel('z [m]');
colorbar;
axis image;
view(2);

subplot(2,2,3);
trisurf(tri,pts(:,1),pts(:,2),pts(:,3),abs(Ez_sct(:,ITX)),...
    'LineStyle','None');
title('E_z^{sct}');
xlabel('x [m]');ylabel('y [m]');zlabel('z [m]');
colorbar;
axis image;
view(2);

subplot(2,2,4);
trisurf(tri,pts(:,1),pts(:,2),pts(:,3),abs(Ez_tot(:,ITX)),...
    'LineStyle','None');
title('E_z^{tot}');
xlabel('x [m]');ylabel('y [m]');zlabel('z [m]');
colorbar;
axis image;
view(2);

figure(2);
probe_theta = atan2(probes(:,2),probes(:,1))*180/pi;
hold off
plot3(probes(:,1),probes(:,2),abs(Ez_inc_p(:,ITX))); hold on;
plot3(probes(:,1),probes(:,2),abs(Ez_sct_p(:,ITX)));
plot3(probes(:,1),probes(:,2),abs(Ez_tot_p(:,ITX)));
grid on;
title(['Field Measurements for tx ' num2str(ITX)]);
legend('Incident','Scattered','Total');

radii = sqrt(sum(centroids.^2,2));
msk = abs(radii-(0.275*LAMBDA)) < (LAMBDA);
t_centroids = centroids(msk,:);
t_theta = atan2(t_centroids(:,2),t_centroids(:,1))*180/pi;
[t_theta,I] = sort(t_theta);
t_centroids = t_centroids(I,:);
t_Ez_tot = Ez_tot(msk,ITX);
t_Ez_tot = t_Ez_tot(I);
figure(3);
ss = scatter(t_theta,abs(t_Ez_tot),...
    'MarkerEdgeAlpha',0.4);
xlabel('Element Centroid Angular Position [deg]');
ylabel('E_Z^T');


phi = linspace(0,pi,100);
ww = 0*phi;
for ii = 1:length(phi)
    wsum = 0;
    cp = cos(phi(ii));
    sp = sin(phi(ii));
    for jj = 1:length(areas)        
        xj = centroids(jj,1);
        yj = centroids(jj,2);
        EjT = Ez_tot(jj,ITX);
        aj = eq_radii(jj);
        J1 = besselj(1,K_B*aj);
        ee = exp(1.0j*K_B*(xj*cp+yj*sp));
        wsum = wsum + ...
            (eps_r(jj)-1)*...
            EjT*...
            aj*...
            J1*...
            ee;
    end
    ww(ii) = pi^2*K_B*abs(wsum)^2;
end

figure(4);
plot(phi*180/pi,ww/LAMBDA);



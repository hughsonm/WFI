clearvars;

ITX = 2;
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
eps_r = k2/(K_B)^2;

probes = dlmread('probe_xyz.txt');
Ez_inc_p = dlmread('Ez_inc_d.txt');
Ez_sct_p = dlmread('Ez_sct_d.txt');
Ez_tot_p = dlmread('Ez_tot_d.txt');


tags = tri(:,end);
tri = tri(:,1:end-1);
tri = tri+1; % Add one because this was made by a zero-based indexing code

f_mesh = figure(1);
subplot(2,2,1);
trisurf(tri,pts(:,1),pts(:,2),pts(:,3),real(eps_r));
title('eps_r');
xlabel('x [m]');ylabel('y [m]');zlabel('z [m]');
colorbar;axis image; view(2);


subplot(2,2,2);
trisurf(tri,pts(:,1),pts(:,2),pts(:,3),abs(Ez_inc(:,ITX)),...
    'LineStyle','None');
title('|E_z^{inc}|');
xlabel('x [m]');ylabel('y [m]');zlabel('z [m]');
colorbar;
axis image;
view(2);

subplot(2,2,3);
trisurf(tri,pts(:,1),pts(:,2),pts(:,3),abs(Ez_sct(:,ITX)),...
    'LineStyle','None');
title('|E_z^{sct}|');
xlabel('x [m]');ylabel('y [m]');zlabel('z [m]');
colorbar;
axis image;
view(2);

subplot(2,2,4);
trisurf(tri,pts(:,1),pts(:,2),pts(:,3),abs(Ez_tot(:,ITX)),...
    'LineStyle','None');
title('|E_z^{tot}|');
xlabel('x [m]');ylabel('y [m]');zlabel('z [m]');
colorbar;
axis image;
view(2);

f_probe = figure(2);
probe_theta = atan2(probes(:,2),probes(:,1))*180/pi;
hold off
plot3(probes(:,1),probes(:,2),abs(Ez_inc_p(:,ITX))); hold on;
plot3(probes(:,1),probes(:,2),abs(Ez_sct_p(:,ITX)));
plot3(probes(:,1),probes(:,2),abs(Ez_tot_p(:,ITX)));
grid on;
title(['Field Measurements for tx ' num2str(ITX)]);
xlabel x;ylabel y;
zlabel('|E|');
legend('Incident','Scattered','Total');


radii = sqrt(sum(centroids.^2,2));
msk = abs(radii-(0.275*LAMBDA)) < (LAMBDA);
t_centroids = centroids(msk,:);
t_theta = atan2(t_centroids(:,2),t_centroids(:,1))*180/pi;
[t_theta,I] = sort(t_theta);
t_centroids = t_centroids(I,:);
t_Ez_tot = Ez_tot(msk,ITX);
t_Ez_tot = t_Ez_tot(I);

f_ez = figure(3);
plot_start_idx = max(1,find(0<t_theta,1)-1);
ss = plot(...
    t_theta(plot_start_idx:end),...
    abs(t_Ez_tot(plot_start_idx:end)),...
    'b*-');%,...
    %'MarkerEdgeAlpha',0.7);
title('Electric Field Distribution Inside Shell');
xlabel('Element Centroid Angular Position [deg]');
ylabel('E_Z^T');


Eimagsq = 1;
if(ITX == 2)
    Eimagsq = abs(besselh(0,2,pi))^2;
end
phi = linspace(0,pi,1000);
echo_width = 0*phi;
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
    echo_width(ii) = pi*pi*K_B*abs(wsum)^2/Eimagsq;
end

f_echo = figure(4);
plot(phi*180/pi,echo_width/LAMBDA);
title('Far-Field Behaviour');
xlabel('\phi [deg]');
ylabel('Echo Width / \lambda');


for ii = 1:size(tri,1)
    tx = pts(tri(ii,:),1);
    ty = pts(tri(ii,:),2);
    ta = polyarea(tx,ty);
    tc = [mean(tx),mean(ty),0];
    assert(abs(ta-areas(ii)) < 1e-8)
    assert(norm(centroids(ii,1:3)-tc) < 1e-7);
end

fig_dest = "Report/figs";
saveas(f_echo,fig_dest + "/Echo" + ITX + ".png");
saveas(f_ez,fig_dest + "/EzInTarget" + ITX + ".png");
saveas(f_mesh,fig_dest + "/Mesh" + ITX + ".png");
saveas(f_probe,fig_dest + "/Probe" + ITX + ".png");




function tightenAxes()
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
end





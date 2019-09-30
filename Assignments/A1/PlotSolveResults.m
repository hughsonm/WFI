clearvars;

ITX = 4;

tri = dlmread('Test_tri.txt');
pts = dlmread('Test_pts.txt');

Ez_inc = dlmread('Ez_inc.txt');
Ez_sct = dlmread('Ez_sct.txt');
Ez_tot = dlmread('Ez_tot.txt');

k2 = dlmread('k2_fgd.txt');

probes = dlmread('probe_xyz.txt');
Ez_inc_p = dlmread('Ez_inc_d.txt');
Ez_sct_p = dlmread('Ez_sct_d.txt');
Ez_tot_p = dlmread('Ez_tot_d.txt');


tags = tri(:,end);
tri = tri(:,1:end-1);
tri = tri+1; % Add one because this was made by a zero-based indexing code

figure(1)
subplot(2,2,1);
trisurf(tri,pts(:,1),pts(:,2),pts(:,3),real(k2));
title('k^2');
xlabel('x [m]');ylabel('y [m]');zlabel('z [m]');
colorbar;axis image; view(2);

subplot(2,2,2);
trisurf(tri,pts(:,1),pts(:,2),pts(:,3),real(Ez_inc(:,ITX)));
title('E_z^{inc}');
xlabel('x [m]');ylabel('y [m]');zlabel('z [m]');
colorbar;
axis image;
view(2);

subplot(2,2,3);
trisurf(tri,pts(:,1),pts(:,2),pts(:,3),real(Ez_sct(:,ITX)));
title('E_z^{sct}');
xlabel('x [m]');ylabel('y [m]');zlabel('z [m]');
colorbar;
axis image;
view(2);

subplot(2,2,4);
trisurf(tri,pts(:,1),pts(:,2),pts(:,3),real(Ez_tot(:,ITX)));
title('E_z^{tot}');
xlabel('x [m]');ylabel('y [m]');zlabel('z [m]');
colorbar;
axis image;
view(2);

figure(2);
hold off
plot3(probes(:,1),probes(:,2),real(Ez_inc_p(:,ITX))); hold on;
plot3(probes(:,1),probes(:,2),real(Ez_sct_p(:,ITX)));
plot3(probes(:,1),probes(:,2),real(Ez_tot_p(:,ITX)));
title(['Field Measurements for tx ' num2str(ITX)]);
legend('Incident','Scattered','Total');



function [fwd_fig,probe_fig] = PlotForwardSolve(FWD_DIR)

if(FWD_DIR(end) ~= '/' || FWD_DIR(end) ~= '\')
    FWD_DIR(end+1) = '/';
end

fwd_eps_r =     ReadCppVectorFromFile([FWD_DIR 'eps_r.txt']);

fwd_tri =       ReadCppMatrixFromFile([FWD_DIR 'tri_tri.txt']);
fwd_pts =       ReadCppMatrixFromFile([FWD_DIR 'tri_pts.txt']);

fwd_tags = fwd_tri(:,end);
fwd_tri = fwd_tri(:,1:end-1);
fwd_tri = fwd_tri+1;

Ez_inc =    ReadCppMatrixFromFile([FWD_DIR 'Ez_inc.txt']);
Ez_sct =    ReadCppMatrixFromFile([FWD_DIR 'Ez_sct.txt']);
Ez_tot =    ReadCppMatrixFromFile([FWD_DIR 'Ez_tot.txt']);

Ez_inc_p =  ReadCppMatrixFromFile([FWD_DIR 'Ez_inc_d.txt']);
Ez_sct_p =  ReadCppMatrixFromFile([FWD_DIR 'Ez_sct_d.txt']);
Ez_tot_p =  ReadCppMatrixFromFile([FWD_DIR 'Ez_tot_d.txt']);

xyz_p =     ReadCppMatrixFromFile([FWD_DIR 'probe_xyz.txt']);

fwd_fig = figure();

ITX = 1;

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

end


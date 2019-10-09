clearvars;
close all;

ITX = 1;
FREQ = 1E9;
CNAUGHT = 2.9979E8;
K_B = (2*pi*FREQ)/CNAUGHT;

LAMBDA = CNAUGHT/FREQ;

tri =       ReadCppMatrixFromFile('Tri.txt');
pts =       ReadCppMatrixFromFile('Points.txt');
alphas =    ReadCppMatrixFromFile('Alphas.txt');
X =         ReadCppVectorFromFile('Chi.txt');
Ez_opt =    ReadCppVectorFromFile('Ez_tot_opt.txt');
% centroids = ReadCppMatrixFromFile('tri_centroids.txt');
% areas =     ReadCppVectorFromFile('tri_areas.txt');
% eq_radii = sqrt(areas/pi);

% Ez_inc =    ReadCppMatrixFromFile('Ez_inc.txt');
% Ez_sct =    ReadCppMatrixFromFile('Ez_sct.txt');
% Ez_tot =    ReadCppMatrixFromFile('Ez_tot.txt');

% k2 =        ReadCppVectorFromFile('k2_fgd.txt');
% eps_r = k2/(K_B)^2;

% probes =    ReadCppMatrixFromFile('ProbePositions.txt');
% Ez_inc_p =  ReadCppMatrixFromFile('Ez_inc_d.txt');
% Ez_sct_p =  ReadCppMatrixFromFile('Ez_sct_d.txt');
% Ez_tot_p =  ReadCppMatrixFromFile('Ez_tot_d.txt');


tags = tri(:,end);
tri = tri(:,1:end-1);
tri = tri+1; % Add one because this was made by a zero-based indexing code


figure();
subplot(3,1,1);
itx = 1;
trisurf(tri,pts(:,1),pts(:,2),0*pts(:,1),real(alphas(:,itx)));
title('Contrast Source');

subplot(3,1,2);
trisurf(tri,pts(:,1),pts(:,2),0*pts(:,1),real(X));

subplot(3,1,3);
trisurf(tri,pts(:,1),pts(:,2),0*pts(:,1),real(Ez_opt(:,itx)));







close all;
clearvars;

tri = dlmread('Test_tri.txt');
pts = dlmread('Test_pts.txt');

Ez_inc = dlmread('Ez_inc.txt');
Ez_sct = dlmread('Ez_sct.txt');

k2 = dlmread('k2_fgd.txt');

tags = tri(:,end);
tri = tri(:,1:end-1);
tri = tri+1; % Add one because this was made by a zero-based indexing code

figure();
subplot(3,1,1);
trisurf(tri,pts(:,1),pts(:,2),pts(:,3),real(k2));
colorbar;axis image; view(2);

subplot(3,1,2);
trisurf(tri,pts(:,1),pts(:,2),pts(:,3),real(Ez_inc(:,1)));
colorbar;
axis image;
view(2);

subplot(3,1,3);
trisurf(tri,pts(:,1),pts(:,2),pts(:,3),real(Ez_sct(:,1)));
colorbar;
axis image;
view(2);


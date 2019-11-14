close all;
clearvars;

FWD_DIR = 'FwdData_A3';
ANTENNA_FILENAME = 'Antennas4A3.txt';
PROBE_FILENAME = 'ProbePositions4A3.txt';
FREQ = 1e9;
OMEGA = 2*pi*FREQ;
C_0 = 2.9979e8;
WAVELENGTH = C_0/FREQ;
k_b = OMEGA/C_0;
k2_b = k_b^2;

fid = fopen(ANTENNA_FILENAME,'r');
tline = fgetl(fid);
sline = split(tline);
assert(strcmpi(sline{1},'x'));
assert(strcmpi(sline{2},'y'));
assert(strcmpi(sline{3},'z'));
assert(strcmpi(sline{4},'magnitude'));
assert(strcmpi(sline{5},'phase'));
assert(strcmpi(sline{6},'style'));
PlaneWaves = zeros(0,3);
while 1
    tline = fgetl(fid);
    if ~ischar(tline), break, end
    sline = split(tline);
    assert(strcmpi(sline{6},'planewave'));
    xy = str2double(sline(1:2)).';
    assert(abs(str2double(sline{3}))<1e-10);
    mp = str2double(sline(4:5)).';
    ri = mp(1)*exp(1.0j*mp(2));
    PlaneWaves = [PlaneWaves;[xy,ri]];
end
fclose(fid);

fid = fopen(PROBE_FILENAME);
tline = fgetl(fid);
sline = split(tline);
assert(strcmpi(sline{1},'x'));
assert(strcmpi(sline{2},'y'));
assert(strcmpi(sline{3},'z'));
Probes = zeros(0,2);
while 1
    tline = fgetl(fid);
    if ~ischar(tline), break, end
    sline = split(tline);
    assert(abs(str2double(sline{3}))<1e-10);
    xy = str2double(sline(1:2)).';
    Probes = [Probes;xy];
end
fclose(fid);

Ez_sct_d = ReadCppMatrixFromFile([FWD_DIR '/Ez_sct_d.txt']);
Ez_inc_d = ReadCppMatrixFromFile([FWD_DIR '/Ez_inc_d.txt']);
Ez_tot_d = ReadCppMatrixFromFile([FWD_DIR '/Ez_tot_d.txt']);

[nrx,ntx] = size(Ez_sct_d);
assert(ntx == size(PlaneWaves,1));
assert(nrx == size(Probes,1));

probe_fig = figure();
for tt=1:ntx
    subplot(3,1,1);
    plot3(Probes(:,1),Probes(:,2),abs(Ez_sct_d(:,tt)));hold on;
    subplot(3,1,2);
    plot3(Probes(:,1),Probes(:,2),abs(Ez_inc_d(:,tt)));hold on;
    subplot(3,1,3);
    plot3(Probes(:,1),Probes(:,2),abs(Ez_tot_d(:,tt)));hold on;
end


ChiHat_pts = zeros(0,3);

ChiHat_fig = figure();

for tt = 1:ntx
    aa = PlaneWaves(tt,1:2).';
    aa = aa/norm(aa);
    for kk = 1:nrx
        rk = Probes(kk,:).';
        rk_dist = norm(rk);
        ss = rk/rk_dist;
        qq = aa-ss;
        cc = k2_b*(1-1.0j)*exp(-1.0j*k_b*rk_dist)/4/sqrt(pi*k_b*rk_dist);
        kx = +qq(1)*k_b;
        ky = +qq(2)*k_b;
        u_samp = Ez_sct_d(kk,tt);
%         if(0<kx)
%             if(0<ky)
%                 % Top-left and bottom-right
%                 ChiHat_pts = [ChiHat_pts;[kx,ky,u_samp/cc]];
%                 ChiHat_pts = [ChiHat_pts;[-kx,-ky,conj(u_samp/cc)]];
%             else
%                 % Bottom-left, ignore.
%             end
%         else
%             if(0<=ky)
%                 % Top-right and bottom-left
%                 ChiHat_pts = [ChiHat_pts;[kx,ky,u_samp/cc]];
%                 ChiHat_pts = [ChiHat_pts;[-kx,-ky,conj(u_samp/cc)]];
%             else
%                 % Bottom-right, ignore
%                 
%             end
%         end
        ChiHat_pts = [ChiHat_pts;[kx,ky,u_samp/cc]];
    end
    
end
ChiHat_full_tri = delaunay(ChiHat_pts(:,1:2));
figure(ChiHat_fig);
subplot(2,1,1);
trisurf(ChiHat_full_tri,ChiHat_pts(:,1),ChiHat_pts(:,2),real(ChiHat_pts(:,3)));
colorbar;
xlabel('k_x');
ylabel('k_y');
subplot(2,1,2);
trisurf(ChiHat_full_tri,ChiHat_pts(:,1),ChiHat_pts(:,2),imag(ChiHat_pts(:,3)));
colorbar;
xlabel('k_x');
ylabel('k_y');
drawnow;

%%
K_RADIUS = k_b*3;

ChiHat_rad = sqrt(sum(ChiHat_pts(:,1:2).^2,2));

ChiHat_pts_restricted = ChiHat_pts(ChiHat_rad<K_RADIUS,:);

tri_hat = delaunay(ChiHat_pts_restricted(:,1:2));

DOM_SIZE = 4*WAVELENGTH;
DOM_NP = 50;

[dom_x,dom_y] = meshgrid(...
    linspace(-DOM_SIZE/2,DOM_SIZE/2,DOM_NP),...
    linspace(-DOM_SIZE/2,DOM_SIZE/2,DOM_NP));
dom_x = dom_x(:);
dom_y = dom_y(:);

tri_chi = delaunay(dom_x,dom_y);
dom_x = dom_x + DOM_SIZE/DOM_NP/4*(rand(length(dom_x),1)-0.5);
dom_y = dom_y + DOM_SIZE/DOM_NP/4*(rand(length(dom_y),1)-0.5);

dom_chi(:) = TIFT(...
    tri_hat,...
    ChiHat_pts_restricted(:,1:2),...
    ChiHat_pts_restricted(:,3),...
    [dom_x,dom_y]);


figure();
subplot(2,1,1);
trisurf(tri_chi,dom_x,dom_y,real(dom_chi));
view(2);
colorbar;
axis image;
subplot(2,1,2);
trisurf(tri_chi,dom_x,dom_y,imag(dom_chi));
view(2);
colorbar;
axis image;

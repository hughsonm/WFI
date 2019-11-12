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

[nrx,ntx] = size(Ez_sct_d);
assert(ntx == size(PlaneWaves,1));
assert(nrx == size(Probes,1));

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
        kx = -qq(1)*k_b;
        ky = -qq(2)*k_b;
        u_samp = Ez_sct_d(kk,tt);
        ChiHat_pts = [ChiHat_pts;[kx,ky,u_samp/cc]];
    end
end
figure(ChiHat_fig);
hold off;
scatter3(ChiHat_pts(:,1),ChiHat_pts(:,2),real(ChiHat_pts(:,3)));
hold on;
scatter3(ChiHat_pts(:,1),ChiHat_pts(:,2),imag(ChiHat_pts(:,3)));
hold off;
xlabel('k_x');
ylabel('k_y');
legend('\Real[\chi]','\Imag[\chi]');
drawnow;


ChiHat_interpolant = scatteredInterpolant(...
    ChiHat_pts(:,1),...
    ChiHat_pts(:,2),...
    ChiHat_pts(:,3),...
    'natural',...
    'none');
%%
RESOLUTION = WAVELENGTH/2;
XY_MAX = WAVELENGTH*2;

RESOLUTION_LIMIT = pi/sqrt(2)/k_b;
if(RESOLUTION < RESOLUTION_LIMIT)
    warning('Requested resolution will lead to numerical abnormalities');
    input('<Enter> to continue\n','s');
end

KXY_BOUND = pi/RESOLUTION;
dxy = pi/KXY_BOUND;
KXY_NPTS = 1+ceil(2*XY_MAX/dxy);



kx_vec = linspace(-KXY_BOUND,KXY_BOUND,KXY_NPTS);
ky_vec = linspace(-KXY_BOUND,KXY_BOUND,KXY_NPTS);

[KX,KY] = meshgrid(kx_vec,ky_vec);
ChiHat = zeros(KXY_NPTS,KXY_NPTS);

ChiHat(:) = ChiHat_interpolant(KX(:),KY(:));
ChiHat(isnan(ChiHat)) = 0;
figure();
subplot(2,1,1)
imagesc([min(kx_vec),max(kx_vec)],...
    [min(ky_vec),max(ky_vec)],...
    real(ChiHat));
colorbar;
subplot(2,1,2);
imagesc([min(kx_vec),max(kx_vec)],...
    [min(ky_vec),max(ky_vec)],...
    imag(ChiHat));colorbar;

Chi = ifft(ifftshift(ChiHat));

fprintf('Spatial Resolution:\t%e lambda\n',dxy/WAVELENGTH);
xy_bound = abs(KXY_NPTS-1)*dxy/2;


figure();
subplot(2,1,1);
imagesc([-xy_bound,xy_bound]/WAVELENGTH,...
    [-xy_bound,xy_bound]/WAVELENGTH,...
    real(Chi));
colorbar;
xlabel('x/\lambda');
ylabel('y/\lambda');
subplot(2,1,2);
imagesc([-xy_bound,xy_bound]/WAVELENGTH,...
    [-xy_bound,xy_bound]/WAVELENGTH,...
    imag(Chi));
colorbar;
xlabel('x/\lambda');
ylabel('y/\lambda');
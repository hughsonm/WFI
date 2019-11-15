close all;
clearvars;

MIN_FREQ = 20e3;
MAX_FREQ = 150e3;

C_0 = 1500;
RHO_0 = 1000;
KAPPA_0 = 1/C_0/C_0/RHO_0;

MIN_WAVELENGTH = C_0/MAX_FREQ;
MAX_WAVELENGTH = C_0/MIN_FREQ;

MIN_K_RAD = 2*2*pi*MIN_FREQ/C_0;
MAX_K_RAD = 2*2*pi*MAX_FREQ/C_0;

NOISE_PCT = 10;

txp = dlmread('A3P4Data/txPositions.dat');

u_scat = load('A3P4Data/uSc.mat');
u_scat = u_scat.uSc;

avg_u_scat = mean(mean(abs(u_scat)));
noise_mags = rand(size(u_scat))*avg_u_scat*NOISE_PCT/100;
noise_angles = rand(size(u_scat))*2*pi;
noise = noise_mags.*exp(1j*noise_angles);
u_scat = u_scat + noise;

ntx = size(txp,1);
nfreq = size(u_scat,2);
assert(size(u_scat,1)==ntx);
frequencies = linspace(MIN_FREQ,MAX_FREQ,nfreq);

ChiHatPts = zeros(ntx*nfreq,3);
ChiHatPtr = 1;

for kk = 1:ntx
    rk = txp(kk,:);
    rk_norm = norm(rk);
    ss = rk/rk_norm;
    for ff = 1:nfreq
        freq = frequencies(ff);
        k_f = 2*pi*freq/C_0;
        cc = -1j*k_f*exp(-1j*2*k_f*rk_norm)/8/pi/rk_norm;
        kx = -2*k_f*ss(1);
        ky = -2*k_f*ss(2);
        u_sample = u_scat(kk,ff);
        ChiHatPts(ChiHatPtr,:) = [kx,ky,u_sample/cc];
        ChiHatPtr = ChiHatPtr+1;
    end
end

ChiHatTri = delaunay(ChiHatPts(:,1:2));
ChiHatTriX = 0*ChiHatTri;ChiHatTriX(:) = ChiHatPts(ChiHatTri(:),1);
ChiHatTriY = 0*ChiHatTri;ChiHatTriY(:) = ChiHatPts(ChiHatTri(:),2);
ChiHatTriRad = sqrt(sum([mean(ChiHatTriX,2),mean(ChiHatTriY,2)].^2,2));
ChiHatTriMsk = MIN_K_RAD < ChiHatTriRad & ChiHatTriRad < MAX_K_RAD/1; 
ChiHatTri = ChiHatTri(ChiHatTriMsk,:);

DOM_SIZE = 28e-2;
DOM_NP = 200;

[dom_x,dom_y] = meshgrid(...
    linspace(-DOM_SIZE/2,DOM_SIZE/2,DOM_NP),...
    linspace(-DOM_SIZE/2,DOM_SIZE/2,DOM_NP));
dom_x = dom_x(:);
dom_y = dom_y(:);

dom_tri = delaunay(dom_x,dom_y);
dom_x = dom_x + DOM_SIZE/DOM_NP/4*(rand(length(dom_x),1)-0.5);
dom_y = dom_y + DOM_SIZE/DOM_NP/4*(rand(length(dom_y),1)-0.5);

dom_chi = TIFT(ChiHatTri,ChiHatPts(:,1:2),ChiHatPts(:,3),[dom_x,dom_y]);
dom_cc = C_0./(sqrt(1+dom_chi));

P4_fig = figure();
subplot(2,2,1);
trisurf(ChiHatTri,ChiHatPts(:,1),ChiHatPts(:,2),real(ChiHatPts(:,3)),...
    'LineStyle','None');
view(2);colorbar;axis equal;
xlabel('k_x');ylabel('k_y');title('Real of F.T. of \chi');
subplot(2,2,2);
trisurf(ChiHatTri,ChiHatPts(:,1),ChiHatPts(:,2),imag(ChiHatPts(:,3)),...
    'LineStyle','None');
view(2);colorbar;axis equal;
xlabel('k_x');ylabel('k_y');title('Imag of F.T. of \chi');
subplot(2,2,3);
trisurf(dom_tri,dom_x,dom_y,real(dom_cc),...
    'LineStyle','None');
view(2);colorbar;axis equal;
xlabel('k_x');ylabel('k_y');title('Real Sound Speed');
subplot(2,2,4);
trisurf(dom_tri,dom_x,dom_y,imag(dom_cc),...
    'LineStyle','None');
view(2);colorbar;axis equal;
xlabel('k_x');ylabel('k_y');title('Imag Sound Speed');

saveas(P4_fig,'Report/figs/P4.png');
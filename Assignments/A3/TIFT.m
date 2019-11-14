function [xx] = TIFT(tri,pp,xhat,rr)
%TIFT Triangular Inverse Fourier Transform
%   Assumes xhat is constant over triangles.
if(size(pp,2) ~= 2)
    if(size(pp,1) ~= 2)
        error('pp should be Nx2');
    else
        pp = pp.';
    end
end

if(size(rr,2) ~= 2)
    if(size(rr,1) ~= 2)
        error('rr should be Nx2');
    else
        rr = rr.';
    end
end
if(size(xhat,2) ~= 1)
    if(size(xhat,1) ~= 1)
        error('xhat should be Nx1');
    else
        xhat = xhat.';
    end
end

xx = zeros(1,size(rr,1));

% tri = delaunay(pp(:,1),pp(:,2));
for ii = 1:size(tri,1)
    tt = tri(ii,:);
    qq = pp(tt(1),:);
    uu = pp(tt(2),:) - qq;
    vv = pp(tt(3),:) - qq;
    J = abs(vv(1)*uu(2)-vv(2)*uu(1));
    xhi = mean(xhat(tt));
    
    
    qdr = qq*(rr.');
    udr = uu*(rr.');
    vdr = vv*(rr.');
    umvdr = (uu-vv)*(rr.');    
    
    x_comp = J*xhi*exp(1j*qdr)/1j./vdr.*(...
        (exp(1j*udr)-exp(1j*vdr))./(1j*umvdr)- ...
        (exp(1j*udr)-1)./(1j*udr));
    
    xx = xx+x_comp;
end
xx = xx.'/((2*pi)^2);




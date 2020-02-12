function f_handle = ViewInversionResults(varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

nargs = length(varargin);

if(nargs<1)
    error('Must provide a data directory as argument, figure handle optional');
elseif(nargs<2)
    data_dir = varargin{1};
    % No figure handle provided. Make one.
    f_handle = figure();
elseif(nargs<3)
    data_dir = varargin{1};
    f_handle = varargin{2};
end

if(~(data_dir(end)=='/' || data_dir(end) == '\'))
    data_dir(end+1) = '/';
end

figure(f_handle);
colormap jet;
max_iter_idx = CountOccurrences(data_dir,'chi_iter_')-1;

pts = ReadCppMatrixFromFile([data_dir 'pts.txt']);
tri = ReadCppMatrixFromFile([data_dir 'tri.txt'])+1; tri = tri(:,1:3);
t2p = BuildTri2PointMap(tri);

for ii = 0:max_iter_idx
    chi_filename = [data_dir 'chi_iter_' num2str(ii) '.txt'];
    Ez_dom_filename = [data_dir 'Ez_tot_iter_' num2str(ii) 'freq_0tx_0.txt'];
    chi = ReadCppVectorFromFile(chi_filename);
    Ez_dom = ReadCppVectorFromFile(Ez_dom_filename);
    subplot(3,2,1);
    trisurf(tri,...
        pts(:,1),pts(:,2),t2p*real(chi),...
        real(chi),...
        'LineStyle','None');
    colorbar;
    xlabel('x');ylabel('y');zlabel('Real(\chi)');
    view(2);
    subplot(3,2,2);
    trisurf(tri,...
        pts(:,1),pts(:,2),t2p*imag(chi),...
        imag(chi),...
        'LineStyle','None');
    colorbar;
    xlabel('x');ylabel('y');zlabel('Imag(\chi)');
    view(2);
    
    subplot(3,2,3);
    trisurf(tri,...
        pts(:,1),pts(:,2),t2p*real(Ez_dom),...
        real(Ez_dom),...
        'LineStyle','None');
    colorbar;
    xlabel('x');ylabel('y');zlabel('Real(E_z^{DOM})');
    view(2);
    
    subplot(3,2,4);
    trisurf(tri,...
        pts(:,1),pts(:,2),t2p*imag(Ez_dom),...
        imag(Ez_dom),...
        'LineStyle','None');
    colorbar;
    xlabel('x');ylabel('y');zlabel('Imag(E_z^{DOM})');
    view(2);
    
    if(ii==0)
        Fs_filename = [data_dir 'Fs.txt'];
        Fs = ReadCppVectorFromFile(Fs_filename);
        subplot(3,2,5:6);
        loglog(Fs);
        xlabel('Iteration');
        ylabel('F_S');
    end
    
    drawnow;
end

end

function n_files = CountOccurrences(dirstring,prefix)
    dirstruct = dir(dirstring);
    n_files = 0;
    for ii = 1:length(dirstruct)
        idir = dirstruct(ii);        
        if(strncmpi(prefix,idir.name,length(prefix)))
            n_files = n_files+1;
        end
    end                
end

function [U,d] = GaussianElimination(A,b)
    M = [A,b];
    nrows = size(A,1);
    ncols = size(A,2);
    assert(nrows==ncols);
    for irow = 1:(nrows-1)
        pp = M(irow,irow);        
        qq = M(irow+1:end,irow);        
        M(irow+1:end,:) = M(irow+1:end,:) - M(irow,:).*qq/pp;                
    end
    U = M(:,1:end-1);
    d = M(:,end);
end

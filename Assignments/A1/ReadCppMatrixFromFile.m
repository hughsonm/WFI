function M = ReadCppMatrixFromFile(filename)
fid = fopen(filename);
dimline = strsplit(fgetl(fid),'\t');
nrows = str2double(dimline{1});
ncols = str2double(dimline{2});
M = zeros(nrows,ncols);
for ii = 1:nrows
    ll = fgetl(fid);
    ss = strsplit(ll,'\t');
    for jj = 1:ncols
        xs = ss{jj};
        if(xs(1) == '(')
            xps = strsplit(strip(strip(xs,'('),')'),',');
            xn = str2double(xps{1}) + 1.0j*str2double(xps{2});
        else
            xn = str2double(xs);
        end
        M(ii,jj) = xn;
    end
end
fclose(fid);
end
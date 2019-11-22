function V = ReadCppVectorFromFile(filename)
fid = fopen(filename);
dimline = strsplit(fgetl(fid),'\t');
nrows = str2double(dimline{1});
V = zeros(nrows,1);
for ii = 1:nrows
    ll = fgetl(fid);
    if(ll(1) == '(')
        lps = strsplit(strip(strip(ll,'('),')'),',');
        nn = str2double(lps{1}) + 1.0j*str2double(lps{2});
    else
        nn = str2double(ll);
    end
    V(ii) = nn;
end
fclose(fid);
end

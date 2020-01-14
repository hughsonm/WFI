function WriteCppMatrixToFile(filename,M)
fid = fopen(filename,'w');

size_str = "";
for dd = size(M)
    if(dd~=1)
        size_str = size_str + dd + sprintf('\t');
    end
end
size_str = size_str{1}(1:end-1);
fprintf(fid,'%s\n',size_str);

all_real = isreal(M);

for ii = 1:size(M,1)
    row_str = "";
    for jj = 1:size(M,2)
        out_val = M(ii,jj);
        if(isnan(out_val)),out_val = 0.0;end
        if(isinf(out_val)),out_val = 0.0;end
        if(all_real)
            row_str = row_str + out_val;
        else
            row_str = row_str + "(" + real(out_val) + "," + imag(out_val) + ")";
        end
        if jj ~= size(M,2)
            row_str = row_str + sprintf('\t');
        end
    end
    row_str = row_str + newline;
    fprintf(fid,'%s',row_str);
end

fclose(fid);
end
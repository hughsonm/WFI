function [Map] = BuildTri2PointMap(tri)
n_tri = size(tri,1);
M_rows = zeros(numel(tri),1);
M_cols = M_rows;
M_vals = M_rows;
M_ptr = 1;

for ii = 1:n_tri
    for jj = tri(ii,:)
        M_rows(M_ptr) = jj;
        M_cols(M_ptr) = ii;
        M_vals(M_ptr) = 1;
        M_ptr = M_ptr+1;
    end
end
Map = sparse(M_rows,M_cols,M_vals);
Map = Map./sum(Map,2);
end


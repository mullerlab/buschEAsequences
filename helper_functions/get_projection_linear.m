function [rho, D, V, prj] = get_projection_linear(spikephases, ndim)

% INPUT
% spikephases = neurons x trials matrix of spike times 
% ndim = dimensions for projection 
% OUTPUT
% rho = correlation matrix
% D = eigenvalues sorted in descending order of magnitude
% V = corresponding eigenvectors (in order) 
% prj = projection onto ndim dimensions

N = size( spikephases, 2); 

rho = ones(N); 
for ii = 1: N
    for jj = 1: ii-1
        x = spikephases(:,ii); y = spikephases(:,jj); 
        xbar = mean( x, 'omitnan'); ybar = mean(y, 'omitnan');
        num = sum( (x-xbar).*(y-ybar) , 'omitnan') ; 
        denom = sqrt( sum( (x-xbar).^2 , 'omitnan' )*sum( (y-ybar).^2 , 'omitnan' ) ); 
        rho(ii,jj) = num/denom; 
    end
end
rho = rho + rho' - 1; 
rho( isnan(rho)) = 0; 

[V,D] = eig(abs(rho)); [D, ind] = sort(abs(diag(D)), 'descend') ; V = V(:,ind); 
Q = V(:,ndim); 
prj = rho*Q; 


end
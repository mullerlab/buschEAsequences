function [r,NUM] = corr_omitnan(x,y, type)

% INPUTS
% x,y = data vectors
% type = 'Pearson' or 'Spearman'

% OUTPUTS
% rho = correlation coefficient 
% NUM = number of entries that are not NAN in both (number of cells being
% compared) 

if length(x) ~= length(y) 
    r = NaN; NUM = NaN; 
else
    n = length(x); 
    
    if strcmp(type, 'Pearson') 
    
    x_mean = mean(x, 'omitnan'); 
    y_mean = mean(y, 'omitnan'); 
    
    num = (x - x_mean).*(y-y_mean); NUM = length(nonzeros(~isnan(num))); 
    num = sum(num, 'omitnan'); 
    
    denom_1 = sum( (x - x_mean).^2 , 'omitnan'); 
    denom_2 = sum( (y - y_mean).^2 , 'omitnan'); 
    denom = sqrt( denom_1 * denom_2 ) ; 
    
    r = num/denom; 
    
    elseif strcmp(type, 'Spearman') 
    
    [~, x_rank] = sort(x, 'descend', 'MissingPlacement', 'last') ;  
    [~, y_rank] = sort(y, 'descend', 'MissingPlacement', 'last') ; 
    
    rank_diff = sum( (x_rank - y_rank).^2 );  
    
    r = 1 - 6*(rank_diff)/(n^3 -n) ; 
    end 
end

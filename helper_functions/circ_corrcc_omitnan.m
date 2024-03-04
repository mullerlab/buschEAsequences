function [rho, pval] = circ_corrcc_omitnan(alpha1, alpha2)
%
% [rho pval] = circ_corrcc(alpha1, alpha2)
%   Circular correlation coefficient for two circular random variables.
%
%   Input:
%     alpha1	sample of angles in radians
%     alpha2	sample of angles in radians
%
%   Output:
%     rho     correlation coefficient
%     pval    p-value
%
% References:
%   Topics in circular statistics, S.R. Jammalamadaka et al., p. 176
%
% PHB 6/7/2008
%

%%%%% ADAPTED FROM THE FOLLOWING FUNCTION TO IGNORE NAN VALUES: 
% Circular Statistics Toolbox for Matlab
% By Philipp Berens, 2009
% berens@tuebingen.mpg.de - www.kyb.mpg.de/~berens/circStat.html


if size(alpha1,2) > size(alpha1,1)
	alpha1 = alpha1';
end

if size(alpha2,2) > size(alpha2,1)
	alpha2 = alpha2';
end

if length(alpha1)~=length(alpha2)
  error('Input dimensions do not match.')
end

% compute mean directions
n = length(alpha1);
alpha1_bar = circ_mean_omitnan(alpha1);
alpha2_bar = circ_mean_omitnan(alpha2);

% compute correlation coeffcient from p. 176
num = sum(sin(alpha1 - alpha1_bar) .* sin(alpha2 - alpha2_bar), 'omitnan');
den = sqrt(sum(sin(alpha1 - alpha1_bar).^2, 'omitnan') .* sum(sin(alpha2 - alpha2_bar).^2, 'omitnan'));
rho = num / den;	

% compute pvalue
l20 = mean(sin(alpha1 - alpha1_bar).^2, 'omitnan');
l02 = mean(sin(alpha2 - alpha2_bar).^2, 'omitnan');
l22 = mean((sin(alpha1 - alpha1_bar).^2) .* (sin(alpha2 - alpha2_bar).^2), 'omitnan');

ts = sqrt((n * l20 * l02)/l22) * rho;
pval = 2 * (1 - normcdf(abs(ts)));
end

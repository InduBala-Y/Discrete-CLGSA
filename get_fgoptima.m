% ******************************************************************************
% * Version: 1.0
% * Last modified on: 21 January, 2016 
% 
% * ****************************************************************************

function [fit] = get_fgoptima(nfunc)
fgoptima = [200.0 1.0 1.0 200.0 1.03163 186.731 1.0 2709.0935 1.0 -2.0 zeros(1,10)];
fit = fgoptima(nfunc);

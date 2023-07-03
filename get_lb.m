% ******************************************************************************
% * Version: 1.0

% * ****************************************************************************

function lb = get_lb(fno)
if (fno == 1 || fno== 2 || fno== 3)
	lb = 0;
elseif (fno== 4)
	lb = -6*ones(1,2);
elseif (fno== 5)
	lb = [-1.9 -1.1];
elseif (fno== 6 || fno== 8)
	lb = -10*ones(1,2);
elseif (fno== 7 || fno== 9)
	lb = 0.25*ones(1,2);
elseif (fno== 10)
	lb = zeros(1,2);
elseif (fno== 11 || fno== 12 || fno== 13)
	dim = 2;
	lb = -5*ones(1,dim);
elseif (fno== 14 || fno== 15)
	dim = 3;
	lb = -5*ones(1,dim);
elseif (fno== 16 || fno== 17)
	dim = 5;
	lb = -5*ones(1,dim);
elseif (fno== 18 || fno== 19)
	dim = 10;
	lb = -5*ones(1,dim);
elseif (fno== 20 )
	dim = 20;
	lb = -5*ones(1,dim);
else
	lb = [];
end

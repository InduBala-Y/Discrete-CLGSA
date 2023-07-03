%
%
% This function gives boundaries and dimension of search space for test functions.
function [down,up,dim]=test_functions_range(nvars)

%If lower bounds of dimensions are the same, then 'down' is a value.
%Otherwise, 'down' is a vector that shows the lower bound of each dimension.
%This is also true for upper bounds of dimensions.

%Insert your own boundaries with a new F_index.

dim=50;
 % Rosenbrock’s Function
    down=-nvars;
    up=nvars;
end

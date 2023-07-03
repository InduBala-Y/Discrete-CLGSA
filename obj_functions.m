
% This function calculates the value of the objective function.
function fit=obj_functions(L,F_index,dim)
max_it=600;
%Insert your own objective function with a new F_index.
F_index=10;
if F_index==10 % Rosenbrock’s Function
    global P_length
  cost=-20*exp(-.2*sqrt(sum(L.^2)/dim))-exp(sum(cos(2*pi.*L))/dim)+20+exp(1)*10*P_length;
    fit=sum(cost);

end

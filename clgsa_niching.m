% ******************************************************************************
% * ****************************************************************************

% Demonstration file on how to use the benchmark suite of the Competition 
% Bellow you can find examples on how to: 
%       i) evaluate a solution, and (niching_func)
%      ii) calculate the number of global optima in a set of solutions (count_goptima)
% 
clear all;

Dims = [1 1 1 2 2 2 2 3 3 2 2 2 2 3 3 5 5 10 10 20]; % dimensionality of benchmark functions
Max_FEs =60000;% [50000*ones(1,5) 200000 200000 400000 400000 200000*ones(1,4) 400000*ones(1,7)]; 
%DO NOT FORGET
noptima = [2 5 1 4 2 18 36 81 216 12 6 8 6 6 8 6 8 6 6 8];
global initial_flag; % the global flag used in test suite 

% for func_num = 1:1
% 	% Set the lower and upper bound for each function
% 	% DO NOT FORGET
% 	initial_flag = 0; % should set the flag to 0 for each run, each function 
% 
% 	% Dimension of the problem
% 	D = Dims(func_num);
% 
% 	% Potential solution 
% 	x = ones(1,D);
% 
% 	% Evaluate the solution
% 	val = niching_func(x, func_num); % fitness evaluation
% 	fprintf('f_%d : f(1...1) = %f\n',func_num,val);
% end
fprintf('---------------------------------------------------------------\n');
fgoptima = [200.0 1.0 1.0 200.0 1.03163 186.731 1.0 2709.0935 1.0 -2.0 zeros(1,10)];
ps=200;
for func_num = 2:2
	initial_flag = 0; % should set the flag to 0 for each run, each function 

	% Dimension of the problem
	D = Dims(func_num);

	% Potential solution 
	x = ones(1,D);

	% Evaluate the solution
	val = niching_func(x, func_num); % fitness evaluation
	fprintf('f_%d : f(1...1) = %f\n',func_num,val);
    runs=1;
    for i=1:runs
	initial_flag = 0; % should set the flag to 0 for each run, each function 
lb=get_lb(func_num);ub=get_ub(func_num);
D=length(lb);acc=[];pos=[];
	% Randomize pulation within optimization bounds 
	% (here dummy initialization within [0,1] ONLY for DEMO)
   p = lb(1) + (ub(1)-lb(1))*rand(ps,D);
 
      varrange=[]; mvden=.1;
  mv=[];
  for i=1:D
      varrange=[varrange;lb(1) ub(1)];
      mv=[mv;(varrange(i,2)-varrange(i,1))/mvden];
  end
if length(mv)==1
 velmaskmin = -mv.*ones(ps, D);     % min vel, psXD matrix
 velmaskmax = mv.*ones(ps, D);
elseif length(mv)==D
    velmaskmin = repmat(forcerow(-mv),ps,1); % min vel
 velmaskmax = repmat(forcerow( mv),ps,1); % max vel
else
 error('Max vel must be either a scalar or same length as prob dimension D');
end

      vel=(repmat(lb,ps,1)-p)+(repmat(ub,ps,1)-p).*rand();%normmat(rand(),[forcecol(xMin-pos(i,j)),xMax-pos(i,j)],1);

      %vel=zeros(ps, D);
     % val =niching_func(p(1,:), func_num); % fitness evaluation
      FEs = ps;
      for i=1:ps
           val(i) =niching_func(p(i,:), func_num); % fitness 
      end
      pbest=p;
      fitPbest=val;
      gen=1; totalGen=Max_FEs; 
      
      min_flag=2;% min_flag=1 for minimization and 2 for maximization
     
      
      while (FEs <= Max_FEs-ps)
          
          
          Pc=0.9;X_b=[];
          dim=D;
for d=1:dim
    if rand()<Pc
        X1=p(ceil(rand*ps),:);
        X2=p(ceil(rand*ps),:);
        F1=niching_func(X1,func_num);% evaluateF(X1,F_index)
        F2=niching_func(X2,func_num); %evaluateF(X2,F_index);
        
        if F1>F2
            X_b = [X_b X1(d)];
        else
            X_b=[X_b X2(d)];
        end 
        
    else 
%          [Fb idz1]=max(fitPbest);
%         X_b=[X_b X(idz1,d)];
    end
end 
                                         % R_X_B=[R_X_B; (X_b)];
                                         % volXb=[volXb; abs(prod(X_b))];
 C_F=niching_func(X_b, func_num);%evaluateF(X_b,F_index)


   
    
          
          
          
          
         % random search, you should add your own method
         % rand('seed', sum(100 * clock));
         LU=[];
          %vel= (repmat(min_var,ps,1)-p)+(repmat(max_var,ps,1)-p).*rand();
        % pbest=p;
         for i=1:ps
            
            X=pbest(i,:);
            fitX= fitPbest(i);
             
        CR = unifrnd(0.4,2.0); F=unifrnd(0.5,1);

                % Implement DE to generate the trial vector
         U = mutation(pbest, [lb; ub], i, ps, D, X, F, CR);
         fitU=niching_func(U, func_num);
          if fitU > fitX
                        pbest(i,:) = U;
                        fitPbest(i) = fitU;
          end
           if C_F>fitU
              pbest(i,:) = X_b;
              fitPbest(i) = C_F;
           end 
          
   LU=[LU U];
         end
         if C_F>fitU
             fitU=C_F;
             U=X_b;
         end 
 G=Gconstant(gen,totalGen); Rnorm=2;Rpower=2;ElitistCheck=1;
[mass]=massCalculation(fitPbest,min_flag);
acceleration=Gfield(mass,p,G,Rnorm,Rpower,ElitistCheck,gen,totalGen);
w=0.9;c1=0.8;c2=1.5;
acc=[acc acceleration(15,1)];
pos=[pos p(15,1)];


                        % Update the velocity of each particle
    vel =w*vel+ c1*(rand(ps, D)).*acceleration+ c2*(rand(ps, D)).* (pbest - p);

                     temp_s = p;

%                         Update the position of each particle
                        p = p + vel;
   vel = ( (vel <= velmaskmin).*velmaskmin ) + ( (vel > velmaskmin).*vel );
   vel = ( (vel >= velmaskmax).*velmaskmax ) + ( (vel < velmaskmax).*vel );   
             for i=1:ps
                        tag1 = find(p(i,:) < lb(1));
                        tag2 = find(p(i,:) > ub(1));

                        if length(tag1)>0
                            p(i,tag1) = (temp_s(1,tag1) + lb(1,tag1)) .* 0.5;
                        end

                        if length(tag2)>0
                            p(i,tag2) = (temp_s(1,tag2) + ub(1,tag2)) .* 0.5;
                        end

             end
             for i=1:ps
           newval(i) =niching_func(p(i,:), func_num); % fitness 
             end
             index=find (fitPbest> newval);
                        p(index,:) = temp_s(index,:);
                        newval(index) = fitPbest(index);
                        fitPbest=newval;
                      
         FEs = FEs + ps;
         gen=gen+1;
	fprintf('f_%d : f(1...1) = %f\n',func_num,fitPbest);
    
      end

	% How many global optima have been found?
	accuracy = 0.00001;
	[count, goptima_found] = count_goptima(pbest, func_num, accuracy);
    aa=[count count/noptima(func_num)];
	fprintf('f_%d, In the current population there are %d global optima!\n', func_num,count);
save('RES_100_partilces.txt','aa','-ascii','-append');
	% Print some stuff :-)
	if count ~=0
		goptima_found;
		for i=1:size(goptima_found,1)
			val = niching_func(goptima_found(i,:), func_num);
			fprintf('F_p: %f, F_g:%f, diff: %f\n', val, fgoptima(func_num), abs(fitPbest - fgoptima(func_num)))
			fprintf('F_p - F_g <= %f : %d\n', accuracy, abs(fitPbest - fgoptima(func_num))<accuracy )
           
		end
    end
    end
end

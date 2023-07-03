function [Xmin,Fmin]=WCA(objective_function,constraints,LB,UB,nvars)
%%
%inputs:
%objective_function is FUNCTION 
%nvars are number variables
%UB upper bound
%LB lower bound
%outputs:
%Xmin is optimum point 
%Fmin is minimum function
%%
%----------------------input----------------------------

%-------------------------------------------------------
Nmax=600;                %number of itarations
nvar=nvars;              %nvars are number variables
N=50;                    %number of water
Nsr=8;                  %number of river
Nwat=N-Nsr;
dmax=1e-5;
%----------------initial river and sea------------------

% AA=[0.100, 0.347, 0.440, 0.539, 0.954, 1.081, 1.174, 1.333, 1.488,...
% 1.764, 2.142, 2.697, 2.800, 3.131, 3.565, 3.813, 4.805, 5.952, 6.572,...
% 7.192, 8.525, 9.300, 10.850, 13.330, 14.290, 17.170, 19.180, 23.680,...
% 28.080, 33.700];
% for iii=22:25

AA=[304.8 406.4 508 609.6 762 1016];
AA=1:length(AA);

for i=1:N
    x(i,:)=LB+((UB-LB).*rand(1,nvar));
%     x(i,:)=LB+((UB-LB).*rand(1));
for j=1:nvar
        dx=abs(AA-x(i,j));
        [mm,index]=min(dx);
        x(i,j)=AA(index);
end


    [c, ceq] = constraints(x(i,:));
    Max_C(i)=max(c);
end

[CS,index]=sort(Max_C,'ascend');
Negativ_C=find(Max_C<0);

if length(Negativ_C)~=0 && length(Negativ_C)>=Nsr+1  
   
    for i=1:length(Negativ_C) 
        xx(i,:)=objective_function(x(Negativ_C(i),:));
    end    
    
    [FS,index1]=sort(xx,'ascend');    
    sr=x(index1(1:Nsr+1),:);  
    
elseif length(Negativ_C)~=0 && length(Negativ_C)<Nsr+1
    
    
    for i=1:length(Negativ_C) 
        xx(i,:)=objective_function(x(Negativ_C(i),:));
    end    
    
    [FS,index1]=sort(xx,'ascend');
    
    sr=x(index1(1:length(Negativ_C)),:);    
    sr(length(Negativ_C)+1:Nsr+1,:)=x(index(length(Negativ_C)+1:Nsr+1),:);
    
elseif length(Negativ_C)==0 
    
    sr=x(index(1:Nsr+1),:);
    
end

SR=sr;



 for i=1:Nsr+1
    cost(i)=objective_function(SR(i,:));
 end
    cs=sort(cost,'ascend')';
for i=1:Nsr+1 
    [cc1,r1]=find(cost==cs(i));
    sr_false(i,:)=SR(r1(1),:);
end


%-------------------------------------------------------
sr_false=sr_false(1:Nsr,:);


cs=cs(1:Nsr+1);
CN=cs-max(cs);    %Cn<0

Pn=abs(CN/(sum(CN)));Pn(Nsr+1)=[];

NCn=round(Nwat*Pn);

if sum(NCn)~=Nwat
    NCn(Nsr)=Nwat-sum(NCn(1:Nsr-1));
end

NCn=sort(NCn,'descend');
sr=sr(1:Nsr,:);


zigma_NCn=[0];
for i=1:Nsr
    zigma_NCn=[zigma_NCn sum(NCn(1:i))]; %use for new sea and river
end

%--------------------initial Water---------------------
lower_bound=repmat(LB,Nwat,1);
upper_bound=repmat(UB,Nwat,1);
x=lower_bound+((upper_bound-lower_bound).*rand(Nwat,nvar)); %initial popultion 

for i=1:Nwat
for j=1:nvar
        dx=abs(AA-x(i,j));
        [mm,index]=min(dx);
        x(i,j)=AA(index);
end
end

WATER=x;



%------------------------------------------------------
sea=sr(1,:);
Locate=1;
river=sr(2:Nsr,:);
WATER;
clear Pn CN SR
%-------------------------------------------------------

for ii=1:Nmax
    
%-------------------------------------------------------

%--------------------move to River---------------------
RIVER=[];
for i=1:Nsr
    R=repmat(sr(i,:),NCn(i),1);
    RIVER=[RIVER;R];
end

C1=2;
Vector =RIVER-WATER; 

WATER =WATER +((C1*rand(size(Vector)).* Vector));

WATER=max(WATER,lower_bound);
WATER=min(WATER,upper_bound);

%--------------------Rain--------------------
dmax=dmax-(dmax/Nmax);
% dmax=dmax/exp(ii/Nmax);
for i=1:Nsr
    if i~=Locate && norm(sr(i,:)-sea)<dmax
       
        c=[1];
%         while max(c)>eps
            new_sr=LB+((UB-LB).*rand(1,nvar));
%             new_sr=LB+((UB-LB).*rand(1));
        for j=1:nvar
        dx=abs(AA-new_sr(1,j));
        [mm,index]=min(dx);
        new_sr(1,j)=AA(index);
        end
            
%             [c, ceq] = constraints(new_sr);
%         end
        sr(i,:)=new_sr;
        
        
        
        
        for j=zigma_NCn(i)+1:zigma_NCn(i+1)
            
            new_WATER=LB+((UB-LB).*rand(1,nvar));
            WATER(j,:)=new_WATER;
            
        end

    elseif i==Locate
        
        for j=zigma_NCn(i)+1:zigma_NCn(i+1)
            if norm( WATER(j,:)-sr(i,:))<eps

                new_WATER=sea+sqrt(0.1).*randn(1,nvar);

%                 new_WATER=sea+ (-1+(2*rand(1,nvar))).* min(min(sea-LB,UB-sea));

%                 new_WATER=LB+((UB-LB).*rand(1,nvar));
                WATER(j,:)=new_WATER;
            end
        end
        
    end
end

WATER=max(WATER,lower_bound);
WATER=min(WATER,upper_bound);
x=WATER;
for i=1:Nwat
for j=1:nvar
        dx=abs(AA-x(i,j));
        [mm,index]=min(dx);
        x(i,j)=AA(index);
end
end
WATER=x;
x=[];

%------------move to Sea and other river---------------

new_sr=[];
for i=1:Nsr
    if i~=Locate 
        
        c=inf;
        [c_sr, ceq] = constraints(sr(i,:));
        
        if max(c_sr)>0
            
        while max(c)>max(c_sr)
            new_sr=(C1.*rand(1,nvar).*(sea-sr(i,:)))+sr(i,:);
            new_sr=max(new_sr,LB);
            new_sr=min(new_sr,UB);
            [c, ceq] = constraints(new_sr);                  
        end 
        
        else
            while max(c)>eps
            new_sr=(C1.*rand(1,nvar).*(sea-sr(i,:)))+sr(i,:);
            new_sr=max(new_sr,LB);
            new_sr=min(new_sr,UB);
            [c, ceq] = constraints(new_sr);
            end
        end
% 
%             if(length(c)~=0 && max(c_sr)>eps && max(c_sr)>max(c))
                sr(i,:)=new_sr;
%             elseif (length(c)~=0 && max(c)<eps )
%                 sr(i,:)=new_sr;
%             else
%                 sr(i,:)=sr(i,:);
%             end

    end
end

x=sr;

for i=1:Nsr
for j=1:nvar
        dx=abs(AA-x(i,j));
        [mm,index]=min(dx);
        x(i,j)=AA(index);
end
end

sr=x;
x=[];

%-----------------NEW Sea and River---------------------
new_sr=[];
for i=1:Nsr
    [c_sr, ceq] = constraints(sr(i,:));    
    for j=zigma_NCn(i)+1:zigma_NCn(i+1)
        [c_water, ceq] = constraints(WATER(j,:));
        if max(c_sr)>eps && max(c_sr)>max(c_water)
            new_sr=WATER(j,:);
            WATER(j,:)=sr(i,:);
            sr(i,:)=new_sr;
        elseif max(c_sr)<eps && max(c_water)<eps && objective_function(sr(i,:))>objective_function(WATER(j,:))
            new_sr=WATER(j,:);
            WATER(j,:)=sr(i,:);
            sr(i,:)=new_sr;
        end
    end
end

%-------------------------------------------------------
% for i=1:Nsr
%      [c, ceq] = constraints(sr(i,:));
%      cost=objective_function(sr(i,:));
%      if cost<objective_function(sea)& max(c)<eps
%         sea=sr(i,:);
%         Locate=i;        
%      end
% end
%          
% Fmin=objective_function(sea);
% Xmin=sea;



for i=1:Nsr
    [c_sr, ceq] = constraints(sr(i,:));  
    [c_sea, ceq] = constraints(sea);
    if max(c_sea)>eps && max(c_sea)>max(c_sr)
        new_sea=sr(i,:);
        sr(i,:)=sea;
        sea=new_sea;
        sr(Locate,:)=sea;        
    elseif max(c_sea)<eps && max(c_sr)<eps && objective_function(sea)>objective_function(sr(i,:))
        new_sea=sr(i,:);
        sr(i,:)=sea;
        sea=new_sea;
        sr(Locate,:)=sea; 
    end 
end
NCn(Locate);
Xmin=sea;
Fmin=objective_function(sea);

% for i=1:Nsr
%      cost1(i)=objective_function(sr(i,:));
% end
% [be,index]=min(cost1);
% sea=sr(index,:);


%-------------------------------------------------------

% FE=cost1;
% [Fmin,index]=min(FE);Xmin=sr(index,:);
% Locate=index;
% sea=Xmin;

for i=1:Nwat
    FC(i)=objective_function(WATER(i,:));
end
%------------------plot--------------------------

% subplot(2,2,1:2)
%  
%    title(['Best Function = ',num2str(Fmin)])
%    xlabel('Iterations')
%    ylabel('function value')
%    plot(ii,Fmin,'--r.','LineWidth',2,...
%                 'MarkerEdgeColor','k')
%    hold on
%    
%  subplot(2,2,3)
%    bar(Xmin)
%    xlabel(['Number of variables = ',num2str(nvar)])        
%    ylabel('global minimum')
%    title('Global Minimum')
%     
%    subplot(2,2,4)
 
   [c, ceq] = constraints(Xmin);
%    if length(c)~=0 
%    plot(ii,max(c),'--r.','LineWidth',2,...
%                 'MarkerEdgeColor','k')
%    xlabel('Iterations')        
%    ylabel('Maximum Constrains')
%    title(['maximum constrains = ',num2str(max(c))])
%    hold on
%    end
%  plot3(WATER(:,1),WATER(:,2),FC,'MarkerFaceColor',[1 1 0],'MarkerEdgeColor',[0 0 0],...
%     'MarkerSize',5,...
%     'Marker','o',...
%     'LineStyle','none',...
%     'Color',[0 0 0])
%    xlabel('x(1) ')        
%    ylabel('x(2)')
%   hold on
%   plot3(sr(:,1),sr(:,2),FE,'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 0],...
%     'MarkerSize',15,...
%     'Marker','h',...
%     'LineStyle','none',...
%     'Color',[0 0 0])
% hold off
%    pause(.01)
 
   
  
%-------------------------------------------------

%--------------------STOP-------------------------
% if abs(F1-Fmin)<eps
%     it=it+1;
% else 
%     it=0;
% end
% 
% if it>=20|ii==Nmax
%     break
% end
% F1=Fmin;

% if abs (Fmin-105.7350856673)<1e-7
%     break
% end

%-------------------------------------------------
 fprintf('%d %15.10f %15.10f\n',ii,Fmin,max(c));
end
% [c, ceq] = constraints(Xmin);
% fprintf('%d %d %15.10f %15.10e\n',iii,ii,Fmin,max(c))
% dmax;
% end
  
    
    

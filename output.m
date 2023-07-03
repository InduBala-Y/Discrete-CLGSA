%% Example 1: Add CV PIPE 
d=epanet('twoloop.inp');
disp('Add a CV pipe')
d.plot;
fromNode = d.getNodeNameID{2};
toNode = d.getNodeNameID{6};
index=d.addLinkPipeCV('CVPipe',fromNode,toNode);
d.plot('links','yes');
d.unload;

%% Example 2: Add CV PIPE with Bin Functions
d=epanet('twoloop.inp');
d.Binplot();
[errcode]=d.addBinCVPipe('CV-P1',fromNode,toNode,1000,10,100);
d.plot('links','yes');
d.unload;

%% Example 3: Add CV PIPE with Bin Functions addBinJunction
d=epanet('twoloop.inp');
d.Binplot('nodes','yes');
% addBinJunction
% arguments: newNodeID,X,Y,ToNodeID,newElevation,newBaseDemand,newDemandPattern,newPipeID,
% newLength,newDiameter,newRoughness,Code
% Add Junction + CV pipe
newID='J1';
[x,y]=ginput(1);
ToNodeID='10'; 
newElevation=500; %ft
newBaseDemand=0;
newDemandPattern='1';
newPipeID='CV-P2';
newLength=1000; %ft
newDiameter=10; %in
newRoughness=100;
Code='CVPIPE';%'CVPIPE', 'PIPE', 'PUMP', 'PRV', 'PSV', 'PBV', 'FCV', 'TCV', 'GPV'
errcode=d.addBinJunction(newID,x,y,newElevation,newBaseDemand,newDemandPattern,newPipeID,...
ToNodeID,newLength,newDiameter,newRoughness,Code);
d.plot('links','yes');
d.unload;

%%CLGSA
function [Fbest,Lbest,BestChart,MeanChart]=GSA(F_index,N,max_it,ElitistCheck,min_flag,Rpower)

 
 Rnorm=2; 
[low,up,dim]=test_functions_range(F_index); 
X=initialization(dim,N,up,low); 
BestChart=[];MeanChart=[];
V=zeros(N,dim);
for iteration=1:max_it
    X=space_bound(X,up,low); 
    fitness=evaluateF(X,F_index); 
    if min_flag==1
    [best best_X]=min(fitness); %minimization.
    else
    [best best_X]=max(fitness); %maximization.
    end        
    if iteration==1
       Fbest=best;Lbest=X(best_X,:);
    end
    if min_flag==1
      if best<Fbest  %minimization.
       Fbest=best;Lbest=X(best_X,:);
      end
    else 
      if best>Fbest  %maximization
       Fbest=best;Lbest=X(best_X,:);
      end
    end
BestChart=[BestChart Fbest];
MeanChart=[MeanChart mean(fitness)];
[M]=massCalculation(fitness,min_flag); 
G=Gconstant(iteration,max_it); 
a=Gfield(M,X,G,Rnorm,Rpower,ElitistCheck,iteration,max_it);
[X,V]=move(X,a,V);
end
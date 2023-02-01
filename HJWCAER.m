function [Xmin,Fmin,SUM_Constraints,NFEs,Elapsed_Time,iters]=HJWCAER(objective_function,constraints,LB,UB,nClust,nVars,nPop,Nsr,dmax,max_it)
% Information

% This code is prepared for single objective function (minimization), constrained, and continuous problems.
% Note that in order to obey the copy-right rules, please cite the published paper properly:
% “Water cycle algorithm - a novel metaheuristic optimization method for solving constrained engineering optimization problems”, Computers & Structures, 110-111 (2012) 151-166.
% "Water cycle algorithm with evaporation rate for solving constrained and unconstrained optimization problems Constraint handling approach", Applied Soft Computing, 30, (2015) 58-71.

% The original code from the above is modified to adapt to clustering problems and it is written for the hybridised version of the ERWCA algorithm, called HJ-WCAER. Please cite the following paper too:
% "Data clustering using hybrid water cycle algorithm and a local pattern search method", Advances in Engineering Software, 153 (2021) 102961.
%%%----------------------------------------------------------------------------------------
% INPUTS:
% objective_function:           Objective function which you wish to minimize or maximize
% constraints                       Given constraints
% LB:                                   Lower bound of a problem
% UB:                                  Upper bound of a problem
% nvars:                               Number of design variables
% Npop                                Population size
% Nsr                                   Number of rivers + sea
% dmax                                Evaporation condition constant
% max_it:                             Maximum number of iterations

% OUTPUTS:
% Xmin:                               Optimum solution
% Fmin:                               Cost/fitness of optimum solution
% SUM_Constraints             Summation of constraint violations
% NFEs:                              Number of function evaluations
% Elapsed_Time                  Elapsed time for optimization process
% Default Values for the WCA
format long g
if (nargin <6 || isempty(nPop)), nPop=50; end
if (nargin <7 || isempty(Nsr)), Nsr=8; end
if (nargin <8 || isempty(dmax)), dmax=1e-6; end    % For unconstrained problems dmax=1e-16 and for constrained problems dmax=1e-05
if (nargin <9 || isempty(max_it)), max_it=30; end
% ------------------------------------------------------------------------
% Create initial population and form sea, rivers, and streams
tic
epss=eps;
N_stream=nPop-Nsr;

ind.position=[];
ind.cost=[];
ind.const=[];

pop=repmat(ind,nPop,1);

for i=1:nPop
    pop(i).position=LB+(UB-LB).*rand(nClust,nVars);
    pop(i).cost=objective_function(pop(i).position);
    c = constraints(pop(i).position);
    pop(i).const=sum(c(c>epss));
end

%----------------Sort to selection of river and sea -----------------------
X_Minus=[];
aa=[pop.const];
COST_MINUS=[pop(aa<=epss).cost];
if ~isempty(COST_MINUS)
    X_Minus=pop(aa<=epss);
    [~,INDEX_M]=sort(COST_MINUS);
    X_Minus=X_Minus(INDEX_M);
end

X_PLUS=[];
SUM_C_PLUS=aa(aa>epss);
COST_PLUS=[pop(aa>epss).cost];
if ~isempty(SUM_C_PLUS)
    AD=unique(SUM_C_PLUS);
    if size(AD,2)==nPop
        X_PLUS=pop(aa>epss);
        [~,INDEX_P]=sort(SUM_C_PLUS);
    else
        X_PLUS=pop(aa>epss);
        [OR,INDEX]=sort(SUM_C_PLUS);
        COST_PLUS=COST_PLUS(INDEX);
        
        kk=1;
        N=0;
        
        for m=1:size(AD,2)
            B=length(find(AD(m)==OR));
            [~,IND]=sort(COST_PLUS(kk:N+B));
            INDEX_P(kk:N+B)=IND;
            kk=B+1;
            N=N+B;
        end
    end
    X_PLUS=X_PLUS(INDEX_P);
end

pop=[X_Minus;X_PLUS];
%------------- Forming Sea ------------------------------------------------
sea=pop(1);
%-------------Forming Rivers ----------------------------------------------
river=pop(2:Nsr);
%------------ Forming Streams----------------------------------------------
stream=pop(Nsr+1:end);
%--------- Designate streams to rivers and sea ----------------------------
cs=[sea.cost;[river.cost]';stream(1).cost];

f=0;
if length(unique(cs))~=1
    CN=cs-max(cs);
else
    CN=cs;
    f=1;
end

NS=round(abs(CN/sum(CN))*N_stream);
if f~=1
    NS(end)=[];
end
NS=sort(NS,'descend');
% ------------------------- Modification on NS -----------------------
i=Nsr;
while sum(NS)>N_stream
    if NS(i)>1
        NS(i)=NS(i)-1;
    else
        i=i-1;
    end
end

i=1;
while sum(NS)<N_stream
    NS(i)=NS(i)+1;
end

if find(NS==0)
    index=find(NS==0);
    for i=1:size(index,1)
        while NS(index(i))==0
            NS(index(i))=NS(index(i))+round(NS(i)/6);
            NS(i)=NS(i)-round(NS(i)/6);
        end
    end
end

NS=sort(NS,'descend');
NB=NS(2:end);

%----------- Main Loop of WCA ---------------------------------------------
disp('******************** Water Cycle Algorithm (WCA)********************');
disp('*Iterations  Function Values  Sum_Const ****************************');
disp('********************************************************************');
FF=zeros(max_it,1);
ER=zeros(1,Nsr-1);
%---------------------- Evaporation rate for different rivers -------------
for k=1:Nsr-1
    ER(k)=(sum(NS(k+1))/(Nsr-1)).*rand;
end

for i=1:max_it
    
    %---------- Moving stream to sea---------------------------------------
    for j=1:NS(1)
        stream(j).position=stream(j).position+2.*rand(1).*(sea.position-stream(j).position);
        
        stream(j).position=min(stream(j).position,UB);
        stream(j).position=max(stream(j).position,LB);
        
        
        stream(j).cost=objective_function(stream(j).position);
        c=constraints(stream(j).position);
        stream(j).const=sum(c(c>epss));
        
        if (sea.const<=epss && stream(j).const<=epss && stream(j).cost<sea.cost)||(sea.const>epss && stream(j).const<=epss)
            
            new_sea=stream(j);
            stream(j)=sea;
            sea=new_sea;
            
        elseif sea.const>epss && stream(j).const>epss && sea.const>stream(j).const
            
            new_sea=stream(j);
            stream(j)=sea;
            sea=new_sea;
            
        end
    end
    
    %---------- Moving Streams to rivers-----------------------------------
    for k=1:Nsr-1
        for j=1:NB(k)
            stream(j+sum(NS(1:k))).position=stream(j+sum(NS(1:k))).position+2.*rand(nClust,nVars).*(river(k).position-stream(j+sum(NS(1:k))).position);
            
            stream(j+sum(NS(1:k))).position=min(stream(j+sum(NS(1:k))).position,UB);
            stream(j+sum(NS(1:k))).position=max(stream(j+sum(NS(1:k))).position,LB);
            
            stream(j+sum(NS(1:k))).cost=objective_function(stream(j+sum(NS(1:k))).position);
            c=constraints(stream(j+sum(NS(1:k))).position);
            stream(j+sum(NS(1:k))).const =sum(c(c>epss));
            
            YES=0;
            if (river(k).const<=epss && stream(j+sum(NS(1:k))).const<=epss && stream(j+sum(NS(1:k))).cost<river(k).cost)||(river(k).const>epss && stream(j+sum(NS(1:k))).const<=epss)
                
                new_river=stream(j+sum(NS(1:k)));
                stream(j+sum(NS(1:k)))=river(k);
                river(k)=new_river;
                YES=1;
                
            elseif river(k).const>epss && stream(j+sum(NS(1:k))).const>epss && river(k).const>stream(j+sum(NS(1:k))).const
                
                new_river=stream(j+sum(NS(1:k)));
                stream(j+sum(NS(1:k)))=river(k);
                river(k)=new_river;
                YES=1;
                
            end
            
            if YES==1
                if  (sea.const<=epss && river(k).const<=epss && river(k).cost<sea.cost)||(sea.const>epss && river(k).const<=epss)
                    
                    new_sea=river(k);
                    river(k)=sea;
                    sea=new_sea;

                elseif sea.const>epss && river(k).const>epss && sea.const>river(k).const
                    
                    new_sea=river(k);
                    river(k)=sea;
                    sea=new_sea;
                    
                end
            end
        end
    end
    %---------- Moving rivers to Sea --------------------------------------
    for j=1:Nsr-1
        river(j).position=river(j).position+2.*rand(nClust,nVars).*(sea.position-river(j).position);
        
        river(j).position=min(river(j).position,UB);
        river(j).position=max(river(j).position,LB);
        
        river(j).cost=objective_function(river(j).position);
        c= constraints(river(j).position);
        river(j).const=sum(c(c>epss));
        
        if (sea.const<=epss && river(j).const<=epss && river(j).cost<sea.cost)||(sea.const>epss && river(j).const<=epss)
            
            new_sea=river(j);
            river(j)=sea;
            sea=new_sea;
 
        elseif sea.const>epss && river(j).const>epss && sea.const>river(j).const
            
            new_sea=river(j);
            river(j)=sea;
            sea=new_sea;
            
        end
    end
    %-------------- Evaporation condition and raining process--------------
    % Check the evaporation condition for rivers and sea
    for k=1:Nsr-1
        if ((norm(river(k).position-sea.position)<dmax) || rand<0.1)
            for j=1:NB(k)
                stream(j+sum(NS(1:k))).position=LB+rand(nClust,nVars).*(UB-LB);
            end
        end
    end
    % Check the evaporation condition for streams and sea
    for j=1:NS(1)
        if ((norm(stream(j).position-sea.position)<dmax))
            stream(j).position=sea.position+sqrt(0.1).*randn(nClust,nVars);
        end
    end
    % Check the evaporation condition among rivers
    for k=1:Nsr-1
        if (exp(-i/max_it)<rand) && (NB(k)<ER(k))
            for j=1:NB(k)
                stream(j+sum(NS(1:k))).position=LB+rand(nClust,nVars).*(UB-LB);
            end
        end
    end
    %----------------------------------------------------------------------
    dmax=dmax-(dmax/max_it);
    
    disp(['Iteration: ',num2str(i),'   Fmin= ',num2str(sea.cost,'%4.5f'),'  Sum_Const= ',num2str(sea.const)]);
    FF(i)=sea.cost;
    %--------Hookes and Jeeves taking over, for further exploitation of the local maxima

    if i>15                % Hookes and Jeeves are only allowed to take over, for after 1st 15 iterations
        a=FF(i);             % Current function values
        b=FF(i-5);          % The 5th previous function values
        if abs(b-a)<1e-8 % Tolerancy value for two function values from ERWCA search, can be changed for different problems of big function values
            nvars=nClust*nVars;  % Number of variables for Hookes-Jeeves method
            startpt=sea.position;  % Starting point of search
            rho=0.5;                     %
            epsilon=1e-8;             % Tolerancy value for Hookes-Jeeves method  
            itermax=1000;             % Maximum iteration number for Hookes-Jeeves method 
            [iters,endpt,fnew]=hooke(nvars, startpt, rho, epsilon, itermax, objective_function);
            Xmin=endpt;   % The function value after the search by Hookes-Jeeves method 
            Fmin=fnew;     % step size of Hookes-Jeeves method  (local search)
            FF(i)=objective_function(Xmin); % the best global function values
            break % the exit the local search if no further improvement is achieved by the set tolerancy value for Hookes-Jeeves method  
        end
    end 

end
% Results and Plot
toc;
Elapsed_Time=toc;
% plot(FF,'LineWidth',2);
% xlabel('Number of Iterations');
% ylabel('Function Values');
NFEs=nPop*max_it;
SUM_Constraints=sea.const;
end


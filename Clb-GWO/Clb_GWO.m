% Competitive learning-based Grey Wolf Optimizer 


function [Alpha_score,Alpha_pos,Convergence_curve]=Clb_GWO(SearchAgents_no,Max_iter,lb,ub,dim,fobj)

% initialize alpha, beta, and delta_pos
Alpha_pos=zeros(1,dim);
Alpha_score=inf; %change this to -inf for maximization problems

Beta_pos=zeros(1,dim);
Beta_score=inf; %change this to -inf for maximization problems

Delta_pos=zeros(1,dim);
Delta_score=inf; %change this to -inf for maximization problems

fitness_new=inf;

%Initialize the positions of search agents
Positions=initialization(SearchAgents_no,dim,ub,lb);

%Fitness evaluation of the initial population
for i=1:size(Positions,1)
	FitVal(i) =fobj(Positions(i,:));
end

%Storing the best score and positions from the fitness evaluation of the initial populatio
BestFitScore = FitVal;
BestFitPos = Positions;



%Used for plotting the converence curves
Convergence_curve=zeros(1,Max_iter);

l=0;% Loop counter

% Main loop
while l<Max_iter
    for i=1:size(Positions,1)  
        
       
        
        % Calculate objective function for each search agent
        fitness = FitVal(i);
        % Update Alpha, Beta, and Delta
        if fitness<Alpha_score
            Alpha_score=fitness; % Update alpha
            Alpha_pos=Positions(i,:);
        end
        
        if fitness>Alpha_score && fitness<Beta_score 
            Beta_score=fitness; % Update beta
            Beta_pos=Positions(i,:);
        end
        
        if fitness>Alpha_score && fitness>Beta_score && fitness<Delta_score 
            Delta_score=fitness; % Update delta
            Delta_pos=Positions(i,:);
        end
    end
    
    
    a=2-l*((2)/Max_iter); % a decreases linearly fron 2 to 0
    
    % Update the Position of search agents including omegas
    for i=1:size(Positions,1)
        for j=1:size(Positions,2)     
                       
            r1=rand(); % r1 is a random number in [0,1]
            r2=rand(); % r2 is a random number in [0,1]
            
            A1=2*a*r1-a; % Equation (2.3)
            C1=2*r2; % Equation (2.4)
            
            D_alpha=abs(C1*Alpha_pos(j)-Positions(i,j)); % Equation (2.5)-part 1
            X1=Alpha_pos(j)-A1*D_alpha; % Equation (2.6)-part 1
                       
            r1=rand();
            r2=rand();
            
            A2=2*a*r1-a; % Equation (2.3)
            C2=2*r2; % Equation (2.4)
            
            D_beta=abs(C2*Beta_pos(j)-Positions(i,j)); % Equation (2.5)-part 2
            X2=Beta_pos(j)-A2*D_beta; % Equation (2.6)-part 2       
            
            r1=rand();
            r2=rand(); 
            
            A3=2*a*r1-a; % Equation (2.3)
            C3=2*r2; % Equation (2.4)
            
            D_delta=abs(C3*Delta_pos(j)-Positions(i,j)); % Equation (2.5)-part 3
            X3=Delta_pos(j)-A3*D_delta; % Equation (2.6)-part 3             
            
            GWOPositions(i,j)=(X1+X2+X3)/3;% Equation (2.7)
            
        end
        
        % Corner Bounding
        GWOPositions(i,:)=min(ub, GWOPositions(i,:));
        GWOPositions(i,:)=max(lb, GWOPositions(i,:));


        
	 % Fitness evaluation of positions upated by the standard GWO procedure     
         FitVal_GWO(i)=fobj(GWOPositions(i,:));

          
          
    end % End of the standard GWO procedure
    

    
% Competetive-Learnig procdure    
    for i=1:SearchAgents_no 	%Loop for all the wolves in the wolfpack
        



        Candidates = [1:i-1 i+1:SearchAgents_no];            % Ensuring that the current member is not the partner
        idx = Candidates(randperm(SearchAgents_no-1,7));     % Selection of seven random omega wolves
        
        R1 = GWOPositions(idx(1),:);                       % Assigning randomly selected solution 1
        R2 = GWOPositions(idx(2),:);                       % Assigning randomly selected solution 2
        R3 = GWOPositions(idx(3),:);                       % Assigning randomly selected solution 3
        R4 = GWOPositions(idx(4),:); 			   % Assigning randomly selected solution 4   		
        R5 = GWOPositions(idx(5),:);			   % Assigning randomly selected solution 5	
 	R6 = GWOPositions(idx(6),:); 			   % Assigning randomly selected solution 6  
        R7 = GWOPositions(idx(7),:);  			   % Assigning randomly selected solution 7    
       
        
        
	
          comp=rand(); % Selection of the competetive learning strategy based on the competetion rate (comp)

          if comp>0.75
              ne=rand; %Selection of the second or the third competetive learning strategy alternatively
              	if ne>0.5
                        % Second competetive learning strategy
			Clb_GWOPositions(i,:)= Alpha_pos+a*rand(1,dim).*(R2-R3); % Equation (3.2)
             	else
			% Third competetive learning strategy
			Clb_GWOPositions(i,:)= Alpha_pos+a*rand(1,dim).*(R4-R5)+a*rand(1,dim).*(R6-R7); % Equation (3.3)
       
             	end
           
            
          else
  		% First competetive learning strategy
		Clb_GWOPositions(i,:)=GWOPositions(i,:)+rand(1,dim).*(R3-Alpha_pos)-rand(1,dim).*(GWOPositions(i,:)+(Beta_pos+Delta_pos)); % Equation (3.1)
          end

 	  % Corner Bounding
       	  Clb_GWOPositions(i,:)=min(ub,   Clb_GWOPositions(i,:));
          Clb_GWOPositions(i,:)=max(lb,   Clb_GWOPositions(i,:));
        
        
  	 % Fitness evaluation of positions upated by the Clb-GWO procedure 
         FitVal_Clb_GWO(i)=fobj(Clb_GWOPositions(i,:));

    end % End of Clb-GWO procedure



%     Convergence_curve(l)=Alpha_score;


     
    % Greedy Selection  

    bett = FitVal_GWO < FitVal_Clb_GWO;        % Comparing the fitness of GWO and Clb_GWO and storing the best fitness                   
    bett_coun = repmat(bett',1,dim);
    
    bettFit = bett.* FitVal_GWO + (1-bett) .* FitVal_Clb_GWO;
    bettPos = bett_coun .* GWOPositions + (1-bett_coun) .* Clb_GWOPositions; 	% Storing the best positions
    
    %% Updating the population
    bett = BestFitScore <= bettFit;                           
    bett_coun = repmat( bett',1,dim);
    
    BestFitScore = bett .* BestFitScore + (1-bett) .* bettFit;
    BestFitPos = bett_coun .* BestFitPos + (1-bett_coun) .* bettPos;
    

    % Returning the updated population
    FitVal = BestFitScore;
    Positions = BestFitPos;
    
  
l=l+1; 

end
end





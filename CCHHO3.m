% Developed in MATLAB R2013b
% Source codes demo version 1.0
% _____________________________________________________

% Main paper:
% Harris hawks optimization: Algorithm and applications
% Ali Asghar Heidari, Seyedali Mirjalili, Hossam Faris, Ibrahim Aljarah, Majdi Mafarja, Huiling Chen
% Future Generation Computer Systems,
% DOI: https://doi.org/10.1016/j.future.2019.02.028
% _____________________________________________________

%  Author, inventor and programmer: Ali Asghar Heidari,
%  PhD research intern, Department of Computer Science, School of Computing, National University of Singapore, Singapore
%  Exceptionally Talented Ph. DC funded by Iran's National Elites Foundation (INEF), University of Tehran
%  03-03-2019

%  Researchgate: https://www.researchgate.net/profile/Ali_Asghar_Heidari

%  e-Mail: as_heidari@ut.ac.ir, aliasghar68@gmail.com,
%  e-Mail (Singapore): aliasgha@comp.nus.edu.sg, t0917038@u.nus.edu
% _____________________________________________________
%  Co-author and Advisor: Seyedali Mirjalili
%
%         e-Mail: ali.mirjalili@gmail.com
%                 seyedali.mirjalili@griffithuni.edu.au
%
%       Homepage: http://www.alimirjalili.com
% _____________________________________________________
%  Co-authors: Hossam Faris, Ibrahim Aljarah, Majdi Mafarja, and Hui-Ling Chen

%       Homepage: http://www.evo-ml.com/2019/03/02/hho/
% _____________________________________________________
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Harris's hawk optimizer: In this algorithm, Harris' hawks try to catch the rabbit.
% HHO  with CrissCross Strategy

% T: maximum iterations, N: populatoin size, CNVG: Convergence curve

function [Rabbit_Location,CNVG]=CCHHO(N,MaxFES,lb,ub,dim,fobj)
% disp('CCHHO is now tackling your problem')
% tic
% initialize the location and Energy of the rabbit
Rabbit_Location=zeros(1,dim);
Rabbit_Energy=inf;

%Initialize the locations of Harris' hawks
X=initialization(N,dim,ub,lb);

CNVG=[];
fitness=ones(N,1)*inf;
fitness_x=ones(N,1)*inf;
fitness_mvc=ones(N,1)*inf;
fitness_mhc=ones(N,1)*inf;
fes=0; % FES counter

t=0; % Loop counter
    for i=1:size(X,1)
        % Check boundries
        FU=X(i,:)>ub;FL=X(i,:)<lb;X(i,:)=(X(i,:).*(~(FU+FL)))+ub.*FU+lb.*FL;
        % fitness of locations
        fitness(i)=fobj(X(i,:));
        fes=fes+1;
        % Update the location of Rabbit
        if fitness(i)<Rabbit_Energy
            Rabbit_Energy=fitness(i);
            Rabbit_Location=X(i,:);
        end
    end
while fes<MaxFES
    
    E1=2*(1-(fes/MaxFES)); % factor to show the decreaing energy of rabbit
    % Update the location of Harris' hawks
    for i=1:size(X,1)
        E0=2*rand()-1; %-1<E0<1
        Escaping_Energy=E1*(E0);  % escaping energy of rabbit
        
        if abs(Escaping_Energy)>=1
            %% Exploration:
            % Harris' hawks perch randomly based on 2 strategy:
            
            q=rand();
            rand_Hawk_index = floor(N*rand()+1);
            X_rand = X(rand_Hawk_index, :);
            if q<0.5
                % perch based on other family members
                X(i,:)=X_rand-rand()*abs(X_rand-2*rand()*X(i,:));
            elseif q>=0.5
                % perch on a random tall tree (random site inside group's home range)
                X(i,:)=(Rabbit_Location(1,:)-mean(X))-rand()*((ub-lb)*rand+lb);
            end
            
        elseif abs(Escaping_Energy)<1
            %% Exploitation:
            % Attacking the rabbit using 4 strategies regarding the behavior of the rabbit
            
            %% phase 1: surprise pounce (seven kills)
            % surprise pounce (seven kills): multiple, short rapid dives by different hawks
            
            r=rand(); % probablity of each event
            
            if r>=0.5 && abs(Escaping_Energy)<0.5 % Hard besiege
                X(i,:)=(Rabbit_Location)-Escaping_Energy*abs(Rabbit_Location-X(i,:));
            end
            
            if r>=0.5 && abs(Escaping_Energy)>=0.5  % Soft besiege
                Jump_strength=2*(1-rand()); % random jump strength of the rabbit
                X(i,:)=(Rabbit_Location-X(i,:))-Escaping_Energy*abs(Jump_strength*Rabbit_Location-X(i,:));
            end
            
            %% phase 2: performing team rapid dives (leapfrog movements)
            if r<0.5 && abs(Escaping_Energy)>=0.5, % Soft besiege % rabbit try to escape by many zigzag deceptive motions
                
                Jump_strength=2*(1-rand());
                X1=Rabbit_Location-Escaping_Energy*abs(Jump_strength*Rabbit_Location-X(i,:));
                
                 fes=fes+1;
                if fobj(X1)<fitness(i) % improved move?
                    X(i,:)=X1;
                else % hawks perform levy-based short rapid dives around the rabbit
                    X2=Rabbit_Location-Escaping_Energy*abs(Jump_strength*Rabbit_Location-X(i,:))+rand(1,dim).*Levy(dim);
                     fes=fes+1;
                    if fobj(X2)< fitness(i) % improved move?
                        X(i,:)=X2;
                    end
                end
            end
            
            if r<0.5 && abs(Escaping_Energy)<0.5, % Hard besiege % rabbit try to escape by many zigzag deceptive motions
                % hawks try to decrease their average location with the rabbit
                Jump_strength=2*(1-rand());
                X1=Rabbit_Location-Escaping_Energy*abs(Jump_strength*Rabbit_Location-mean(X));
                
                fes=fes+1;
                if fobj(X1)<fitness(i) % improved move?
                    X(i,:)=X1;
                else % Perform levy-based short rapid dives around the rabbit
                    X2=Rabbit_Location-Escaping_Energy*abs(Jump_strength*Rabbit_Location-mean(X))+rand(1,dim).*Levy(dim);
                     fes=fes+1;
                    if fobj(X2)<fitness(i) % improved move?
                        X(i,:)=X2;
                    end
                end
            end
            FU=X(i,:)>ub;FL=X(i,:)<lb;X(i,:)=(X(i,:).*(~(FU+FL)))+ub.*FU+lb.*FL;
            % Cross
            Bvc = randperm(dim);
            Mvc(i,:) = X(i,:);
            %normalization
            
            Boundary_no= size(ub,2); % numnber of boundaries
            % If the boundaries of all variables are equal and user enter a signle
            % number for both ub and lb
            if Boundary_no==1
                Mvc(i,:) = (Mvc(i,:)-lb)/(ub-lb);
            end
            % If each variable has a different lb and ub
            if Boundary_no>1
                for j=1:dim
                    ub_j=ub(j);
                    lb_j=lb(j);
                    Mvc(i,j)=(Mvc(i,j)-lb_j)/(ub_j-lb_j);
                end
            end
            
            p2 = 0.6;  %p2取0.2到0.8之间
            for j=1:(dim/2)
                p = unifrnd(0,1);  %生成服从均匀分布的0-1的随机数
                if  p<p2
                    no1= Bvc(2*j-1);
                    no2 = Bvc(2*j);
                    r = unifrnd(0,1);  %生成服从均匀分布的0-1的随机数
                    Mvc(i,no1)=r*Mvc(i,no1)+(1-r)*Mvc(i,no2);
                end
            end
            
            Boundary_no= size(ub,2); % numnber of boundaries
            % If the boundaries of all variables are equal and user enter a signle
            % number for both ub and lb
            if Boundary_no==1
                Mvc(i,:) = Mvc(i,:)*(ub-lb)+lb;
            end
            % If each variable has a different lb and ub
            if Boundary_no>1
                for j=1:dim
                    ub_j=ub(j);
                    lb_j=lb(j);
                    Mvc(i,j)=(ub_j-lb_j)*Mvc(i,j)+lb_j;
                end
            end
            % Check boundries
            FU=Mvc(i,:)>ub;FL=Mvc(i,:)<lb;Mvc(i,:)=(Mvc(i,:).*(~(FU+FL)))+ub.*FU+lb.*FL;
            FU=X(i,:)>ub;FL=X(i,:)<lb;X(i,:)=(X(i,:).*(~(FU+FL)))+ub.*FU+lb.*FL;
            % fitness of locations
            fitness_x(i)=fobj(X(i,:));
            fitness_mvc(i)=fobj(Mvc(i,:));
            fes=fes+2;
            %%
        end
    end
        fitness_mix=[fitness_x;fitness_mvc];
        X_mix=[X;Mvc];
        [fitness_mix,index]=sort(fitness_mix, 'ascend');
        for i =1:size(X,1)
            X(i,:)=X_mix(index(i),:);
            fitness(i)=fitness_mix(i);
        end
        if fitness(1)<Rabbit_Energy
            Rabbit_Energy=fitness(1);
            Rabbit_Location=X(1,:);
        end 
    
    % Criss
    Mhc = zeros(size(X,1),dim);
    Bhc = randperm(size(X,1));
    for i=1:(size(X,1)/2)
        %           p = unifrnd(0,1);  %生成服从均匀分布的0-1的随机数
        %           p<p1
        no1= Bhc(2*i-1);
        no2 = Bhc(2*i);
        for j=1:dim
            r1 = unifrnd(0,1);  %生成服从均匀分布的0-1的随机数
            r2 = unifrnd(0,1);  %生成服从均匀分布的0-1的随机数
            c1 = (rand(1)*2)-1; %生成服从均匀分布的-1到1的随机数
            c2 = (rand(1)*2)-1;
            Mhc(no1,j)=r1*X(no1,j)+(1-r1)*X(no2,j)+c1*(X(no1,j)-X(no2,j));
            Mhc(no2,j)=r2*X(no2,j)+(1-r2)*X(no1,j)+c2*(X(no2,j)-X(no1,j));
        end
    end
    
    for i=1:size(X,1)
        % Check boundries
        FU=Mhc(i,:)>ub;FL=Mhc(i,:)<lb;Mhc(i,:)=(Mhc(i,:).*(~(FU+FL)))+ub.*FU+lb.*FL;
        % fitness of locations
        fitness_mhc(i)=fobj(Mhc(i,:));
        fes=fes+1;
    end
    fitness_mix=[fitness;fitness_mhc];
    X_mix=[X;Mhc];
    [fitness_mix,index]=sort(fitness_mix, 'ascend');
    for i =1:size(X,1)
        X(i,:)=X_mix(index(i),:);
        fitness(i)=fitness_mix(i);
    end
    if fitness(1)<Rabbit_Energy
        Rabbit_Energy=fitness(1);
        Rabbit_Location=X(1,:);
    end 
    
    
    
    t=t+1;
    CNVG(t)=Rabbit_Energy;
    %    Print the progress every 100 iterations
    %    if mod(t,100)==0
    %        display(['At iteration ', num2str(t), ' the best fitness is ', num2str(Rabbit_Energy)]);
    %    end
end
% toc
end

% ___________________________________
function o=Levy(d)
beta=1.5;
sigma=(gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);
u=randn(1,d)*sigma;v=randn(1,d);step=u./abs(v).^(1/beta);
o=step;
end


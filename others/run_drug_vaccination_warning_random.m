function [track_exposed, exposed_listSum, who_infect_list, ...
    isolation_Num,  testingNumAll, ...
    hospital_incidence, death_incidence, V2_list, immunityRatio, ...
    vaccine_second_ever_list, treatmentNumAll] = ...
    run_drug_vaccination_warning_random(effi_drug_hosp, EverInfectedPro, VEDs, VEIs, VESs,...
    Ages100, t_Cases, adherence, ...%symptom_blocking,...
    endDays, beta_others, beta_house, house_neighs, neighborsOthers, ...
    sigma, gamma, gammaA, prop_sym, epsilon, ...
    Ages5, baselineFlag, treat_rate, vac_times_ages, trans3, drug_type, eta_Treat_Days, vaccine_ratio_ini, onlyTreatUnvacOlder, iRun)

if exist('onlyTreatUnvacOlder')==0
    onlyTreatUnvacOlder = 0;
end

isolation_Num = 0; testingNumAll=0;
% track_treateddays = [];
immunityRatio = 0;

if treat_rate == -1
    vaccine_ratio_ini = [0, 0, 0, 0, 0];
else
%     vaccine_ratio_ini = 0.7;%0.75;
%     vaccine_ratio_ini = 0;%0.75;
end
% popImmunity = 0.003;
popImmunity = 0.01;
% popImmunity = 0.02;
% popImmunity = 1/3;


% relativeVED1 = 1.-VEDs(1,:);
relativeVED2 = 1.-VEDs(2,:);
% relativeVED3 = 1.-VEDs(3,:);
relativeVEDR = 1.-VEDs(4,:);

% relativeVEI1 = 1.-VEIs(1,:);
relativeVEI2 = 1.-VEIs(2,:);
% relativeVEI3 = 1.-VEIs(3,:);
relativeVEIR = 1.-VEIs(4,:);

% for death
% relativeVES1 = 1.-VESs(1,:);
relativeVES2 = (1-VESs(2,:))./(1-VEDs(2,:));
% relativeVES3 = 1.-VESs(3,:);
relativeVESR = (1.-VESs(4,:))./(1-VEDs(4,:));

% if noEfficacy2Infect==1
%     relativeVEI2 = ones(length(relativeVEI2),1);
%     relativeVEIR = ones(length(relativeVEIR),1);
% end
%% symptomatic case hospitalization ratio
% https://www.nature.com/articles/s41562-020-0931-9#Fig1
% [10] Verity, R. et al. Estimates of the severity of coronavirus disease 2019: a model-based analysis. The Lancet Infectious Diseases 20, 669–677 (2020).
p_sym_hosp = [0, 0.025, 2.672, 9.334, 15.465]*0.01;
eta = 1/5.9;
eta_Treat = 1/eta_Treat_Days;
eta_Treatleft = 1/(5.9-1/eta_Treat);


%% https://wwwnc.cdc.gov/eid/article/26/10/20-1702_article
% death rate on hospitalized individuals, age specific
% H to D: deathRate*death_mu
% H to R: (1-deathRate)*gamma_H
% deathRate = [0.0390, 0.1208, 0.0304, 0.1049, 0.2269];
deathRate = [ 0.0048    0.0048    0.0568    0.1082    0.1615];

death_mu = 0.12821;
gamma_H = 0.0912409;

%% N=network size, T=transmissibility, tau=recovery probability
N = length(Ages5);
vaccine_second_ever_list = zeros(N,1);

randNums = rand(N, endDays, 11);

% exposed_list is a list that contains zeros and ones.
suscept_list = ones(N,1); %Initialize list with zeros indicating all nodes are susceptibles
exposed_list = zeros(N,1); %Initialize list with zeros indicating all nodes are susceptibles
presymp_list = zeros(N,1); %Initialize list with zeros indicating all nodes are susceptibles
sympU_list = zeros(N,1); %Initialize list with zeros indicating none of the nodes are recovered at t=0
sympT_list = zeros(N,1); %Initialize list with zeros indicating none of the nodes are recovered at t=0

% ratio_RiskByTreat = 1-0.891; % https://www.nejm.org/doi/full/10.1056/NEJMoa2118542  0.109
ratio_RiskByTreat = 1-effi_drug_hosp; % https://www.nejm.org/doi/full/10.1056/NEJMoa2118542  0.109

% symp_list
asym_list = zeros(N,1); %Initialize list with zeros indicating none of the nodes are recovered at t=0
recover_list = zeros(N,1); %Initialize list with zeros indicating none of the nodes are recovered at t=0
% V1_list = zeros(N,1); %Initialize list with zeros indicating none of the nodes are recovered at t=0
V2_list = zeros(N,1); %Initialize list with zeros indicating all nodes are susceptibles
symp_ever_list = zeros(N,1);
hospital_list = zeros(N,1); %Initialize list with zeros indicating none of the nodes are recovered at t=0
death_list = zeros(N,1); %Initialize list with zeros indicating none of the nodes are recovered at t=0

if baselineFlag==1
    p_sym2treat = 0*ones(5,1);
else
    p_sym2treat = treat_rate*ones(5,1);
end

%%
hospital_incidence = zeros(N,1);
V1_time_list = 10^5*ones(N,1);
V2_time_list = 10^5*ones(N,1);
who_infect_list = 0;%10^5*ones(N,1);
track_exposed = 10^5*ones(N,1); %tracks time at infection for each node.
track_treated = 10^5*ones(N,1); %tracks time at infection for each node.
track_YT = 10^5*ones(N,1);

% track_hospital = 10^5*ones(N,1); 
% track_death = 10^5*ones(N,1);  

% p_randperm = randperm(N);
% temp = ceil(N*vaccine_ratio_ini);
% p_V2 = p_randperm( (end-temp+1):end);
% vaccine_ratio_ini = [0, 0.406, 0.761, 0.875, 0.096];% [0, 40.6%, 76.1%, 87.5%, 96%]
p_V2 = [];
for i=1:length(Ages5)
    if rand<vaccine_ratio_ini(Ages5(i))
        p_V2 = [p_V2; i];
    end
end
suscept_list(p_V2) = 0;
V2_list(p_V2) = 1;

% V2_time_list(p_V2) = 0;
% vac_times_ages, Ages100
for i=1:length(p_V2)
    temp = vac_times_ages{Ages100(p_V2(i))};
    V2_time_list(p_V2(i)) = temp(randi(length(temp),1)) - datenum('2022-01-29');
end
recover_time_list = 10^5*ones(N,1);


% EverInfectedPro
temp_ps= [];
for i=1:length(EverInfectedPro)
    if i==1;    temp = find(Ages5==1 | Ages5==2); end
    if i>1;     temp = find(Ages5==(i+1)); end
    
    temp_p = temp(randperm( length(temp), ceil(length(temp)*EverInfectedPro(i)) ));
    temp_ps=[temp_ps; temp_p];
end
suscept_list(temp_ps) = 0;
recover_list(temp_ps) = 1;
% recover_time_list(temp_ps) = 0;
recover_time_list(temp_ps) = t_Cases(randperm(length(t_Cases), length(temp_ps))) - datenum('2022-01-29'); % t_Cases


% p_zero = randperm(N, 10);% ceil(N*0.01)); %choose a random node to infect
% temp = find(suscept_list>-1); 
% temp = find(suscept_list>-1); 
temp = find(recover_list==1 | suscept_list==1); 
p_zero = temp(randperm(length(temp), ceil(length(temp)*popImmunity) ));
% p_zero = randperm(N, ceil(N*popImmunity) );
length(p_zero)/N;

track_time = 0;
suscept_list(p_zero) = 0;
exposed_list(p_zero) = 1;
track_exposed(p_zero) = 0;
recover_list(p_zero) = 0;
V2_list(p_zero) = 0;

treatmentNumAll = zeros(N,1);
exposed_listSum = zeros(endDays,1);

ToH = -1*ones(N,1); % if -1: not applied, if 0, to R, if 1 to H
ToHrate = zeros(N,1);

ToH2 = -1*ones(N,1); % if -1: not applied, if 0, to R, if 1 to H
% ToHrate2 = zeros(N,1);

sympT_list_new = zeros(N,1);
sympU_list_new = zeros(N,1);

for ic=1:endDays
    
    vaccine_second_ever_list(find(V2_list==1)) = 1;
%     immunityRatio_org = length(find(recover_list==1) )...% + length(find(Trecover_list==1) ) ...
%         + length(find(V1_list==1) ) + length(find(V2_list==1) );% ...
%     
%     immunityRatio = immunityRatio_org./N;
    
    exposed_listSum(ic) = sum(sympT_list_new) + sum(sympU_list_new);
    sympT_list_new = zeros(N,1);
    sympU_list_new = zeros(N,1);

    track_time = track_time + 1;
    %%
    V1_list_new = zeros(N,1);

    %% # new infected
    % infect the susceptibles with the probability 1-exp(-bk) where k is the number of infeced neighbors.
    exposed_list_new = zeros(N,1);
    NodeList = find(suscept_list==1);
    for i=1:length(NodeList)
        node=NodeList(i);
        pTemp = 0;
        
        %% infected_neighborsSEIR(trans3, track_exposed, track_treated, drug_type, track_time, neighborsAll, node, infected_list)
        pTemp = pTemp+(1- exp(-beta_house *infected_neighborsSEIR(trans3, track_exposed, track_treated, drug_type, track_time,  house_neighs, node,  exposed_list, asym_list, presymp_list, sympU_list, sympT_list, recover_list)));
        pTemp = pTemp+(1- exp(-beta_others*infected_neighborsSEIR(trans3, track_exposed, track_treated, drug_type, track_time,  neighborsOthers, node,  exposed_list, asym_list, presymp_list, sympU_list, sympT_list, recover_list)));
        
        %% mainly for Re estimationw when without testing.
        % So do not consdier the impact of isolated people no matter they
        % have infectiousness or not.
        temp_rand = randNums(node, ic, 1); % randNums = rand(N, endDays, 10);
        if temp_rand<pTemp % && rel_beta_isolated<=0
            exposed_list_new(node)=1;
            temp01 = house_neighs{node};
            temp01 = temp01(find(track_exposed(temp01)<10^5));
            temp = track_treated(temp01)-track_exposed(temp01);
            temp(find(temp>size(trans3,2))) = size(trans3,2);
            
            temp11 = [];
            for temp_i = 1:length(temp01)
                temp_trackT = track_time-track_exposed(temp01(temp_i));
                if temp_trackT > size(trans3, 1)
                    temp_trackT = size(trans3, 1);
                end
                temp11 = [temp11; trans3(temp_trackT, temp(temp_i), drug_type)];
            end
            
            temp02 = neighborsOthers{node};
            temp02 = temp02(find(track_exposed(temp02)<10^5));
            temp = track_treated(temp02)-track_exposed(temp02);
            temp(find(temp>size(trans3,2))) = size(trans3,2);
            
            temp12 = [];
            for temp_i = 1:length(temp02)
                temp_trackT = track_time-track_exposed(temp02(temp_i));
                if temp_trackT > size(trans3, 1)
                    temp_trackT = size(trans3, 1);
                end
                
                temp12 = [temp12; trans3(temp_trackT, temp(temp_i), drug_type)];
            end
            
            temp = [temp01, temp02];
            temp1 = [beta_house.*temp11; beta_others.*temp12];
            
            temp2 = temp(find( mnrnd(1, temp1./sum(temp1)) == 1));
%             who_infect_list(node)= temp2;
            
            if (ismember(temp2, p_zero))
                who_infect_list = who_infect_list+1;
            end
        end
    end
    
    %% Those vaccinated
    for ii=[2, 4]
%         if ii==1; NodeList = find(V1_list==1); end
        if ii==2; NodeList = find(V2_list==1); end
        if ii==4; NodeList = find(recover_list==1); end
        
        for i=1:length(NodeList)
            node=NodeList(i);
            
            temp_beta = 1;
            if ii==2% | ii==4;
                if V2_time_list(node)<10^5
%                     track_time-V2_time_list(node)
                    temp_beta = relativeVEI2(track_time-V2_time_list(node));
                end
            end
            
            if ii==4
                if recover_time_list(node)<10^5
                    temp_beta = relativeVEIR(track_time-recover_time_list(node));
                end
            end
            %% infected_neighborsSEIR(trans3, track_exposed, track_treated, drug_type, track_time, neighborsAll, node, infected_list)
            pTemp = 0;
            pTemp = pTemp+(1- exp(-beta_house *infected_neighborsSEIR(trans3, track_exposed, track_treated, drug_type, track_time,  house_neighs,  node, exposed_list, asym_list, presymp_list, sympU_list, sympT_list, recover_list)));
            pTemp = pTemp+(1- exp(-beta_others*infected_neighborsSEIR(trans3, track_exposed, track_treated, drug_type, track_time,  neighborsOthers,node,exposed_list, asym_list, presymp_list, sympU_list, sympT_list, recover_list)));
            pTemp = pTemp*temp_beta;
            
            temp_rand = randNums(node, ic, 2); % randNums = rand(N, endDays, 10);
            if temp_rand<pTemp;
                exposed_list_new(node)=1;
                
                temp01 = house_neighs{node};
                temp01 = temp01(find(track_exposed(temp01)<10^5));
                temp = track_treated(temp01)-track_exposed(temp01);
                temp(find(temp>size(trans3,2))) = size(trans3,2);
                
                temp11 = [];
                for temp_i = 1:length(temp01)
                    temp_trackT = track_time-track_exposed(temp01(temp_i));
                    if temp_trackT > size(trans3, 1)
                        temp_trackT = size(trans3, 1);
                    end
                    temp11 = [temp11; trans3(temp_trackT, temp(temp_i), drug_type)];
                end
                
                temp02 = neighborsOthers{node};
                temp02 = temp02(find(track_exposed(temp02)<10^5));
                temp = track_treated(temp02)-track_exposed(temp02);
                temp(find(temp>size(trans3,2))) = size(trans3,2);
                
                temp12 = [];
                for temp_i = 1:length(temp02)
                    temp_trackT = track_time-track_exposed(temp02(temp_i));
                    if temp_trackT > size(trans3, 1)
                        temp_trackT = size(trans3, 1);
                    end
                    temp12 = [temp12; trans3(temp_trackT, temp(temp_i), drug_type)];
                end
                
                temp = [temp01, temp02];
                temp1 = [beta_house.*temp11; beta_others.*temp12];
                
                temp2 = temp(find( mnrnd(1, temp1./sum(temp1)) == 1));
%                 who_infect_list(node)= temp2;
                if (ismember(temp2, p_zero))
                    who_infect_list = who_infect_list+1;
                end
                
            end
            
        end
    end
    
    % All the infected nodes recover with a probailty
    asym_list_new = zeros(N,1);
    presymp_list_new = zeros(N,1);
%     symp_list_new = zeros(N,1);
    recovered_list_new = zeros(N,1);
    asymrecovered_list_new = zeros(N,1);
    hospital_list_new = zeros(N,1);
    death_list_new = zeros(N,1);
    
    treated_list_new = zeros(N,1);
    %%
    NodeList = find(exposed_list==1);
    for i=1:length(NodeList)
        node=NodeList(i);
        temp_rand = randNums(node, ic, 3);
        if temp_rand<sigma % epsilon
            temp_prop_sym = prop_sym;
            
            if V2_time_list(node)>=10^5 & V1_time_list(node)<10^5;
%                 temp_prop_sym = (1-0.01*psi1)*prop_sym;
                temp_psi  = relativeVED1(track_time-V1_time_list(node)+1);
                temp_prop_sym = temp_psi*prop_sym;
            end
            
            if V2_time_list(node)<10^5;
%                 temp_prop_sym = (1-0.01*psi2)*prop_sym;
                temp_psi  = relativeVED2(track_time-V2_time_list(node)+1);
                temp_prop_sym = temp_psi*prop_sym;
            end

            if recover_time_list(node)<10^5;
                temp_psi  = relativeVEDR(track_time-recover_time_list(node)+1);
                temp_prop_sym = temp_psi*prop_sym;
            end

            temp_rand = randNums(node, ic, 4);
            if temp_rand < temp_prop_sym
                presymp_list_new(node)=1;
            else
                asym_list_new(node)=1;
            end
        end
    end
    
    %% Those in pre-symptomatic
    % ToH2 = -1*ones(N,1); 
    % ToHrate2 = zeros(N,1);
    NodeList = find(presymp_list==1);
    for i=1:length(NodeList)
        node=NodeList(i);
        if ToH2(node) == -1
            pTreat = p_sym2treat(Ages5(node));
            flag_TreatUnvacOlder = 1;
            if onlyTreatUnvacOlder==1
                if Ages5(node) <=4 & vaccine_second_ever_list(node) == 0
                    flag_TreatUnvacOlder = 0;
                end
            end
            temp_rand = randNums(node, ic, 5);
            if temp_rand<pTreat & flag_TreatUnvacOlder
                ToH2(node) = 1; % treat
            else
                ToH2(node) = 0; % not treat
            end
        else
            temp_rand = randNums(node, ic, 6);
            if temp_rand < epsilon%sigma
                if ToH2(node) == 1
                    sympT_list_new(node)=1;
                    temp = 0;
                    
                    rng(node+ic+iRun);
                    while rand<eta_Treat
                        temp = temp+1;
                    end
                    track_treated(node) = track_time+temp;


                    temp = 0;
                    rng(node+ic+iRun);
                    while rand<eta_Treatleft
                        temp = temp+1;
                    end
                    track_YT(node) = track_treated(node)+temp;

                else
                    sympU_list_new(node)=1;
                end
            end
        end
    end
    
    %% ToH = -1*ones(N,1); % if -1: not applied, if 0, to R, if 1 to H
    % Those in symptomatic
    NodeList = find(sympU_list==1);
    ratio = 1;
    for i=1:length(NodeList)
        node=NodeList(i);
        if ToH(node) == -1
            pYR = 1-p_sym_hosp(Ages5(node))*ratio;
            pYH = p_sym_hosp(Ages5(node))*ratio;
            
            temp_rand = randNums(node, ic, 7);
            randV = temp_rand;
            if randV < pYR
                ToH(node) = 0;
                ToHrate(node) = gamma;
            else
                randV = randV - pYR;
                if randV <= pYH
                    ToH(node) = 1; % hospital
                    ToHrate(node) = eta;
                end
            end
        end
        
        temp_rand = randNums(node, ic, 10);
        if temp_rand<ToHrate(node)
            if ToH(node) == 0 % recover
                recovered_list_new(node)=1;
            end
            if ToH(node) == 1 % hospital
                hospital_list_new(node)=1;
            end
        end
    end

    %
    NodeList = find(sympT_list==1);
    ratio = ratio_RiskByTreat;
    for i=1:length(NodeList)
        node=NodeList(i);
        if track_YT(node)<=track_time
            pYR = 1-p_sym_hosp(Ages5(node))*ratio(Ages5(node));
            temp_rand = randNums(node, ic, 11);
            if temp_rand<pYR
                recovered_list_new(node)=1;
            else
                hospital_list_new(node)=1;
            end
        end
    end
    %%
    % H to D: deathRate*death_mu
    % H to R: (1-deathRate)*gamma_H
    NodeList = find(hospital_list==1);
    for i=1:length(NodeList)
        node=NodeList(i);

        temp_pdeath = 1;
        if V2_time_list(node)<10^5;
            temp_pdeath= relativeVES2(track_time-V2_time_list(node)+1);
        end

        if recover_time_list(node)<10^5;
            temp_pdeath= relativeVESR(track_time-recover_time_list(node)+1);
        end
        
        temp_rand = randNums(node, ic, 8);
        temprand = temp_rand;
        tempRandt = (temp_pdeath*deathRate(Ages5(node))*death_mu ...
                + (1-temp_pdeath*deathRate(Ages5(node)))*gamma_H );


        if temprand< tempRandt
            if temprand< (temp_pdeath*deathRate(Ages5(node))*death_mu)
                death_list_new(node)=1;
            else
                recovered_list_new(node)=1;
            end
        end
    end
    
    %%
    NodeList = find(asym_list==1);
    for i=1:length(NodeList)
        node=NodeList(i);
        temp_rand = randNums(node, ic, 9);
        if temp_rand<gammaA
            asymrecovered_list_new(node)=1;
        end
    end
    
    
    %% Store the time of recovery for the newly nodes in list
    suscept_list((exposed_list_new==1)) = 0;
    suscept_list((V1_list_new==1)) = 0;
    
    %%
    V2_list((exposed_list_new==1)) = 0;
    V2_list((asymrecovered_list_new==1)) = 0;
    V2_list((recovered_list_new==1)) = 0;

    % update the original exposed_list
    exposed_list((exposed_list_new==1)) = 1;
    exposed_list((presymp_list_new==1)) = 0;
%     exposed_list((symp_list_new==1)) = 0;
    exposed_list((sympU_list_new==1)) = 0;
    exposed_list((sympT_list_new==1)) = 0;
    exposed_list((asym_list_new==1)) = 0;
    exposed_list((recovered_list_new==1)) = 0;
    
    % update orignial presymp_list
    presymp_list((exposed_list_new==1)) = 0;
    presymp_list((presymp_list_new==1)) = 1;
%     presymp_list((symp_list_new==1)) = 0;
    presymp_list((sympU_list_new==1)) = 0;
    presymp_list((sympT_list_new==1)) = 0;
    presymp_list((asym_list_new==1)) = 0;
    presymp_list((recovered_list_new==1)) = 0;
    
    % update orignial symp_list
    sympU_list((exposed_list_new==1)) = 0;
    sympU_list((presymp_list_new==1)) = 0;
    sympU_list((sympU_list_new==1)) = 1;
    sympU_list((sympT_list_new==1)) = 0;
    sympU_list((asym_list_new==1)) = 0;
    sympU_list((recovered_list_new==1)) = 0;
    sympU_list((hospital_list_new==1)) = 0;

    % update orignial symp_list
    sympT_list((exposed_list_new==1)) = 0;
    sympT_list((presymp_list_new==1)) = 0;
    sympT_list((sympU_list_new==1)) = 0;
    sympT_list((sympT_list_new==1)) = 1;
    sympT_list((asym_list_new==1)) = 0;
    sympT_list((recovered_list_new==1)) = 0;
    sympT_list((hospital_list_new==1)) = 0;

    %
    asym_list((exposed_list_new==1)) = 0;
    asym_list((presymp_list_new==1)) = 0;
    asym_list((sympT_list_new==1)) = 0;
    asym_list((sympU_list_new==1)) = 0;
    asym_list((asym_list_new==1)) = 1;
    asym_list((asymrecovered_list_new==1)) = 0;
    
    %
    hospital_list((hospital_list_new==1)) = 1;
    hospital_list((death_list_new==1)) = 0;
    hospital_list((recovered_list_new==1)) = 0;

    death_list((death_list_new==1)) = 1;

    recover_list((recovered_list_new==1)) = 1;
    recover_list((asymrecovered_list_new==1)) = 1;
    recover_list((exposed_list_new==1)) = 0;

    %
    hospital_incidence(hospital_list_new==1) = hospital_incidence(hospital_list_new==1)+1;
    treatmentNumAll = treatmentNumAll+ sympT_list_new;

    everTreated(treated_list_new==1) = 1;
    symp_ever_list(find(sympU_list==1)) = 1;
    symp_ever_list(find(sympT_list==1)) = 1;

    %% Store the time of infection for the newly infected nodes in track_exposed list
    track_exposed(find(exposed_list_new==1))= track_time;
    track_treated(find(exposed_list_new==1))= 10^5;
    
    recover_time_list((recovered_list_new==1)) = track_time;
    recover_time_list((asymrecovered_list_new==1)) = track_time;
    V2_time_list((recovered_list_new==1)) = 10^5;
    V2_time_list((asymrecovered_list_new==1)) = 10^5;


%     for i=1:length(sympT_list_new)
%         if sympT_list_new(i)==1 & sympU_list_new(i)==1
%             1;
%         end
%     end
end
death_incidence = death_list;

who_infect_list = who_infect_list/length(p_zero);
  

% track_exposed = 10^5*ones(N,1); %tracks time at infection for each node.
% track_treated = 10^5*ones(N,1); %tracks time at infection for each node.
% track_YT = 10^5*ones(N,1);
% hist(track_treated(find(track_treated<10^5)))
% hist(track_exposed(find(track_treated<10^5)))
% track_treateddays

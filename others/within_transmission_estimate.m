
Triangular_dis_i = 1
temp_str = ['_', num2str(net_house_num_i), num2str(EverInfectedPros_i), ...
    num2str(Rel_Vir_Infect_i), num2str(Triangular_dis_i), '_', num2str(VEI_natural_i)];

virusType = 1;
load trans3_Nirmatrelvir.mat
trans3(:,:,3) = trans3(:,:,2);
if strcmp(Rel_Vir_Infect, 'log10')
end
if strcmp(Rel_Vir_Infect, 'log10-proportional')
    for i=1:size(trans3,1)
        for j=1:size(trans3,2)
            for k=1:size(trans3,3)
                if trans3(i,j,k)<6
                    trans3(i,j,k) = 0;
                else
                    trans3(i,j,k) = trans3(i,j,k)-6;
                end
            end
        end
    end
end

if strcmp(Rel_Vir_Infect, 'threshold')
    for i=1:size(trans3,1)
        for j=1:size(trans3,2)
            for k=1:size(trans3,3)
                if trans3(i,j,k)<6
                    trans3(i,j,k) = 0;
                else
                    trans3(i,j,k) = 10;
                end
            end
        end
    end
end


if strcmp(Rel_Vir_Infect, 'sigmoid')
    temp = temp_sigmoid;
    for i=1:size(trans3,1)
        for j=1:size(trans3,2)
            for k=1:size(trans3,3)
                if trans3(i,j,k)>0
                    [M1, M2] = min(abs(trans3(i,j,k)-temp(:,1) ));
                    trans3(i,j,k) = temp(M2,2)*10;
                end
            end
        end
    end
end


beta_house = 0.35/sum(trans3(:,end,1));
adherence = 1;


drug_type = 2; % Pax
pricePerTest = 5;
pricePerFullVac = 12*2;

load(['US_', num2str(NodeNum), temp_str, '_R0_20221121_vaccination.mat']);
Rtei = randperm(length(Rte_Range));
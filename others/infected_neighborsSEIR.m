function infected_neighborsget = infected_neighborsSEIR(trans3, track_exposed, ...
    track_treated, drug_type, track_time, neighborsAll, node, ...
    exposed_list, asym_list, presymp_list, sympU_list, sympT_list, recover_list)%, infected_list)

neighbors = neighborsAll{node};
infected_neighborsget = 0;
for i=1:length(neighbors)
    if track_exposed(neighbors(i)) < 10^5 && ...
            (exposed_list(neighbors(i))==1 || asym_list(neighbors(i))==1 ...
            || presymp_list(neighbors(i))==1 || sympU_list(neighbors(i))==1 || sympT_list(neighbors(i))==1 || recover_list(neighbors(i))==1 )
        temp_trackT = track_time-track_exposed(neighbors(i));
        if temp_trackT > size(trans3, 1)
            temp_trackT = size(trans3, 1);
        end
        
        if (track_treated(neighbors(i))-track_exposed(neighbors(i)) )>size(trans3, 2)
            infected_neighborsget = infected_neighborsget + trans3(temp_trackT, end, drug_type);
        else
%             if track_treated(neighbors(i))-track_exposed(neighbors(i))<=0
%                 0;
%             end
            infected_neighborsget = infected_neighborsget + trans3(temp_trackT, track_treated(neighbors(i))-track_exposed(neighbors(i)), drug_type);
        end
    end
end


VEI_first = [zeros(1,14), 0.46*ones(1, 1500)];
VED_first = [zeros(1,14), 0.40*ones(1, 1500)];
VES_first = [zeros(1,14), 0.848*ones(1, 1500)];

VEI_second = [0.64*ones(1, 7),  0.64*ones(1,(9-2+1)*7), 0.428*ones(1, (14-10+1)*7),  0.22*ones(1, 1500)];
VED_second = [0.669*ones(1, 7), 0.672*ones(1,(4-2+1)*7), 0.55*ones(1, (9-5+1)*7),    0.457*ones(1, 1500)];
VES_second = [0.919*ones(1, 7), 0.919*ones(1, (9-2+1)*7), 0.919*ones(1, (14-10+1)*7), 0.919*ones(1, 1500)];


VEI_third = [0.928*ones(1, 14), 0.767*ones(1, (9-2+1)*7),0.468*ones(1, (14-10+1)*7), 0.238*ones(1, (19-15+1)*7), 0.100*ones(1, 1500)];
VED_third = [0.974*ones(1, 14), 0.958*ones(1,(9-2+1)*7), 0.843*ones(1, (14-10+1)*7), 0.650*ones(1, (19-15+1)*7), 0.428*ones(1, 1500)];
VES_third = [0.994*ones(1, 14), 0.991*ones(1,(9-2+1)*7), 0.978*ones(1, (14-10+1)*7), 0.951*ones(1, (19-15+1)*7), 0.837*ones(1, 1500)];

VEIs = [VEI_first(1:1500); VEI_second(1:1500); VEI_third(1:1500); VEI_natural(1:1500)];
VEDs = [VED_first(1:1500); VED_second(1:1500); VED_third(1:1500); VED_natural(1:1500)];
VESs = [VES_first(1:1500); VES_second(1:1500); VES_third(1:1500); VES_natural(1:1500)];

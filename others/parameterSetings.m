load vac_times_ages.mat % vac_times_ages
load t_Cases.mat % t_Cases

net_house_nums = [1000, 5000];
Rte_Range = [1.2, 1.5, 1.7, 2, 3];
p_sym_hosp = [0, 0.025, 2.672, 9.334, 15.465]*0.01;
deathRate = [ 0.0048    0.0048    0.0568    0.1082    0.1615];
prop_sym = 0.75;
EverInfectedPros = {143/129*[35.49, 54.86, 45.66, 43.65,  32.36]*0, 143/129*[35.49, 54.86, 45.66, 43.65,  32.36]*0.01};
Rel_Vir_Infects = {'log10'; 'log10-proportional'; 'threshold'; 'sigmoid';};

Triangular_18_49 =  makedist("Triangular","a",0.48,"b", 3*0.59-0.48-0.71,"c", 0.71);
Triangular_50_64 =  makedist("Triangular","a",0.34,"b", 3*0.40-0.34-0.48,"c", 0.48);
Triangular_65    =  makedist("Triangular","a",0.48,"b", 3*0.53-0.48-0.58,"c", 0.58);

vaccine_ratio_ini = [0, 0.424, 0.772, 0.896, 0.974];
onlyTreatUnvacOlder = 0;
treat_rates = [0, 0.2, 0.5 0.8 1];
eta_Treat_Dayset = [3];
pricePerTreatment = 530;
willingnessRange = 100000;
VEI_natural_i = 1;
VEI_natural = [0.831*ones(1, 7), 0.831*ones(1,(9-2+1)*7), 0.731*ones(1, (14-10+1)*7),  0.633*ones(1, 1500)];
VED_natural = [0.921*ones(1, 7), 0.921*ones(1,(4-2+1)*7), 0.892*ones(1, (9-5+1)*7),    0.87*ones(1, 1500)];
VES_natural = [0.981*ones(1, 7), 0.981*ones(1, (9-2+1)*7), 0.981*ones(1, (14-10+1)*7), 0.981*ones(1, 1500)];

endDays=250;
VaccineRates = 0.01*[0.01;0.01;0.01;0.02;0.01;0.02;0.02;0.01;0.01;0.03;0;0.02;0.04;0.04;0.02;0;0.02;0.01;0.02;0.04;0.03;0.03;0;0.03;0.03;0.03;0.03;0.04;0;0.02;0.01;0.01;0.03;0;0.01;0;0.01;0;0.01;0;0;-0.01;0.01;0.02;0.04;0.05;0.05;0.07;0.06;0.01;0.1;0.05;0.05;0.03;0.09;0.07;0;0.16;0.11;0.11;0.1;0.15;0.15;0.03;0.19;0.14;0.14;0.14;0.16;0.14;0.01;0.03;0.12;0.11;0.23;0.26;0.17;0.2;0.23;0.24;0.19;0.24;0.2;0.2;0.18;0.21;0.13;0.18;0.12;0.1;0.1;0.05;0.09;0.11;0.05;0.15;0.06;0.03;0.02;0.03;0.03;0.01;0.02;0;0.02;0;0;0;-0.02;0.01;0;0;0.02;0.01;0.02;0.05;0.04;0.08;0.06;0.07;0.08;0.16;0.06;0.16;0.07;0.18;0.18;0.12;0.23;0.09;0.11;0.14;0.15;0.16;0.32;0.12;0.16;0.14;0.19;0.18;0.2;0.21;0.14;0.18;0.2;0.14;0.16;0.26;0.18;0.25;0.2;0.54;0.36;0.21;0.37;0.35;0.4;0.16;0.13;0.24;0.06;0.21;0.02;0.23;0.1;0.06;0.07;0.07;0.29;0.15;0.03;0.08;0.23;0.05;0.09;0.05;0.09;0.15;0.07;0.05;0.14;0.09;0.14;0.05;0.2;0.09;0.15];
VaccineRates(find(VaccineRates<0)) = 0;
VaccineRates = [VaccineRates; zeros(200,1)];
VaccineRates = VaccineRates(1:endDays);


load networkHousehold_meanUS_one5000_2021.mat
G = networkHousehold.edges;
Ages5 = networkHousehold.ages;
Ages100= networkHousehold.ages100;
Degree = networkHousehold.degrees;
degreeAve = (mean(Degree.^2)-mean(Degree) )/mean(Degree);
NodeNum = length(Degree);
neighbors = networkHousehold.neighs;

neighborsOthers = neighbors;
for i=1:length(neighbors)
    temp = neighborsOthers{i};
    neighborsOthers{i} = temp(~ismember(temp, house_neighs{i}));
end

testingTotDay = endDays;
iRunTot = 2;%100;
YLL_setting;

settingVEs_Paxovir;
hospitalCost = [21847, 21847, 21847*1/32 + 19785*4/32 + 17219*8/32 + 19543*10/32+ 21736*9/32, 21736*1/15 + 24012*10/15 + 21373*4/15, 21373*6/(79-65+1) + 17094*(79-71+1)/(79-65+1)]';
sigma = 1/3;
gammaA = 1/9;
eta_Treat_Days = 3;
gamma = 1/(4);
epsilon = 1/2;

temp_sigmoid = [2.54000000000000,0.00312077200000000;2.60000000000000,0.00312077300000000;2.66000000000000,0.00312075700000000;2.72000000000000,0.00311990700000000;2.78000000000000,0.00312579800000000;2.84000000000000,0.00395780700000000;2.90000000000000,0.00468096500000000;2.96000000000000,0.00468199500000000;3.02000000000000,0.00468116500000000;3.08000000000000,0.00468115300000000;3.14000000000000,0.00468112000000000;3.20000000000000,0.00468325900000000;3.26000000000000,0.00471704700000000;3.32000000000000,0.00641993300000000;3.38000000000000,0.00743759700000000;3.44000000000000,0.00781132900000000;3.50000000000000,0.00780273200000000;3.56000000000000,0.00779909200000000;3.62000000000000,0.00763744200000000;3.68000000000000,0.00932977200000000;3.74000000000000,0.0110709700000000;3.80000000000000,0.0109227840000000;3.86000000000000,0.0109516760000000;3.92000000000000,0.0123399420000000;3.98000000000000,0.0140100110000000;4.04000000000000,0.0140372500000000;4.10000000000000,0.0143999780000000;4.16000000000000,0.0172159160000000;4.22000000000000,0.0170981570000000;4.28000000000000,0.0199999540000000;4.34000000000000,0.0203257110000000;4.40000000000000,0.0224235570000000;4.46000000000000,0.0234411020000000;4.52000000000000,0.0253133410000000;4.58000000000000,0.0265006210000000;4.64000000000000,0.0300392100000000;4.70000000000000,0.0327524660000000;4.76000000000000,0.0324872560000000;4.82000000000000,0.0358883950000000;4.88000000000000,0.0388703520000000;4.94000000000000,0.0422848190000000;5,0.0453392410000000;5.06000000000000,0.0482208700000000;5.12000000000000,0.0521509590000000;5.18000000000000,0.0570107460000000;5.24000000000000,0.0610274630000000;5.30000000000000,0.0651876930000000;5.36000000000000,0.0706960110000000;5.42000000000000,0.0748755380000000;5.48000000000000,0.0812059700000000;5.54000000000000,0.0874397450000000;5.60000000000000,0.0936003550000000;5.66000000000000,0.101274937000000;5.72000000000000,0.108983733000000;5.78000000000000,0.116706694000000;5.84000000000000,0.125651024000000;5.90000000000000,0.135309208000000;5.96000000000000,0.149682862000000;6.02000000000000,0.160554664000000;6.08000000000000,0.167973408000000;6.13993252400000,0.181645515000000;6.20006747600000,0.189529487000000;6.26000000000000,0.203911212000000;6.32000000000000,0.217519497000000;6.38000000000000,0.232832269000000;6.44000000000000,0.248416745000000;6.50000000000000,0.263589468000000;6.56000000000000,0.281285644000000;6.62000000000000,0.299032372000000;6.68000000000000,0.317865306000000;6.74000000000000,0.336038670000000;6.80000000000000,0.355948075000000;6.86000000000000,0.375860103000000;6.92000000000000,0.394834393000000;6.98000000000000,0.416472834000000;7.04000000000000,0.436989799000000;7.10000000000000,0.458661589000000;7.16000000000000,0.480033188000000;7.22000000000000,0.501893142000000;7.28000000000000,0.523158349000000;7.34000000000000,0.543927428000000;7.40000000000000,0.564727953000000;7.46000000000000,0.586383005000000;7.52000000000000,0.606398881000000;7.58000000000000,0.626613466000000;7.64000000000000,0.646449528000000;7.70000000000000,0.665731714000000;7.76000000000000,0.683204982000000;7.82000000000000,0.700824957000000;7.88000000000000,0.719030284000000;7.94000000000000,0.736097608000000;8,0.749622428000000;8.06000000000000,0.766871346000000;8.12000000000000,0.779641257000000;8.18000000000000,0.794492520000000;8.24000000000000,0.805802107000000;8.30000000000000,0.817559631000000;8.36000000000000,0.830792490000000;8.42000000000000,0.842129749000000;8.48000000000000,0.851185010000000;8.54000000000000,0.860821320000000;8.60000000000000,0.869026244000000;8.66000000000000,0.879060511000000;8.72000000000000,0.887162417000000;8.78000000000000,0.893109001000000;8.84000000000000,0.900908245000000;8.90000000000000,0.908723571000000;8.96000000000000,0.913103683000000;9.02000000000000,0.920116520000000;9.08000000000000,0.925425135000000;9.14000000000000,0.930227314000000;9.20000000000000,0.934579641000000;9.26000000000000,0.939169485000000;9.32000000000000,0.943858833000000;9.38000000000000,0.947429223000000;9.44000000000000,0.950207263000000;9.50000000000000,0.953525293000000;9.56000000000000,0.956750418000000;9.62000000000000,0.959973621000000;9.68000000000000,0.961943686000000;9.74000000000000,0.963216932000000;9.80000000000000,0.965692327000000;9.86000000000000,0.969072357000000;9.92000000000000,0.969768685000000];
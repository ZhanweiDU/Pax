% YLL
% Years load and Ages
% https://www.cdc.gov/nchs/data/nvsr/nvsr66/nvsr66_04.pdf
Years_expect = [78.9,0;78.3,1;74.4,5;69.5,10;64.5,15;59.7,20;54.9,25;50.2,30;45.4,35;40.7,40;36.1,45;31.7,50;27.4,55;23.3,60;19.4,65;15.7,70;12.3,75;9.20,80;6.70,85;4.60,90;3.20,95;2.30,100];
% SeasonYear and five age Groups: AGE 0-4	AGE 5-24	AGE 25-49	AGE 50-64	AGE 65
Years_expect5 = [mean(Years_expect(1:2,1)); % 1 for 0-4
    sum([Years_expect(3:4,1)*5; Years_expect(5,1)*3])/13; % 5 to 17: 5,10 and 15 for 5-9, 10-14, 15-19
    sum([Years_expect(5,1)*2; Years_expect(6:11,1)*5])/32; % 18 to 49: 20:5:45 for 20-24, 25-29，30-34， 35-40， 41-44， 45-49
    mean(Years_expect(12:14,1)); % 50 to 64: 50:5:50 for 50-54， 55-59， 60-64
    mean(Years_expect(15:17,1))]; % 65+: :5:75 for 65-69, 70-74, 75+

for i=1:length(Years_expect5)
    temp = round(Years_expect5(i));
    temp1 = 0;
    for j=1:(temp)
        temp1 = temp1+0.97^(j-1);
    end
    Years_expect5(i) = temp1;
end
% Years_expect5

% https://sites.cns.utexas.edu/sites/default/files/cid/files/covid-19_analysis_for_austin_march2020.pdf
% Overall: [0.0016, 0.0049, 0.084, 1.000, 3.371], IFR: infected fatality ratio, age specific (%)
Fatality_ratio = [0.0016, 0.0049, 0.084, 1.000, 3.371]*0.01;

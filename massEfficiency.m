clear all
clc
%% import Data (TABLE 5-6 Fuel Economy Effect of Tare Weight Reduction by Truck Class)
importMassEffData

%% Post Process
truckClasss = numel(massEffData.TruckType);
x = linspace(1,truckClasss,truckClasss);
yBar = zeros(2,truckClasss);
figure1 = figure(1); clf,
for i = 1:truckClasss
    y = [ massEffData.massEffMIn(i) massEffData.massEfmax(i)];
    if y(1) > y(2)
        y = [y(2) y(1)];
    end
    yBar(:,i) = y';
    plot(x(i),y,'ks','MarkerSize',10)
    R=rectangle('Position',[i-0.2 yBar(1,i) 0.4 yBar(2,i)-yBar(1,i)],FaceColor=[0.4660 0.6740 0.1880],EdgeColor='none');
    hold on
end
ylabel('[%]')
xticks(x)
xticklabels(string(massEffData.TruckType))
grid on
ax = gca; 
ax.TickLabelInterpreter = "latex";
ax.FontSize = 20;

pictureWidth = 21; % This one is picture width
hw_ratio = 0.8;  
set(figure1,'Units','centimeters','Position',[60 50 pictureWidth hw_ratio*pictureWidth])
[t,s] = title('paylod-to-total mass ratio $\eta_m$','FontWeight','bold','Interpreter','latex');

%% ghandariz etl. 2021  https://doi.org/10.3390/en14113221
curbWeightGhandariz  =   [10000 33000];  
grossWeightGhandariz =   [25000 80000]; 

payLoadGhandariz = grossWeightGhandariz - curbWeightGhandariz; 

massEffGhandariz = 100 - ((curbWeightGhandariz ./ grossWeightGhandariz).*100); 
hold on  
p = plot(x(end-1:end),massEffGhandariz,"pentagram"); 
p.MarkerFaceColor = [0.9290 0.6940 0.1250]; 
p.MarkerSize = 10;
p.MarkerEdgeColor = [0.9290 0.6940 0.1250]; 

text(x(end-1),massEffGhandariz(1),'\leftarrow ghandariz etl. 2021')
text(x(end),massEffGhandariz(2),'\leftarrowghandariz etl. 2021')

%% M. Fries etl 2017. doi: 10.1109/EVER.2017.7935872
%M. Fries, M. Kruttschnitt and M. Lienkamp, 
% "Multi-objective optimization of a long-haul truck hybrid operational strategy and a predictive powertrain control system,"
% 2017 Twelfth International Conference on Ecological Vehicles and Renewable Energies (EVER), 
% Monte Carlo, Monaco, 2017, pp. 1-7, doi: 10.1109/EVER.2017.7935872

curbWeightFries = 14000; %14 ton
grossWeightFries = curbWeightFries + 25000; 
payLoadFries =  grossWeightFries - curbWeightFries; 
massEffFries = 100 - ((curbWeightFries/grossWeightFries)*100); 

p = plot(x(end),massEffFries,"pentagram"); 
p.MarkerFaceColor = [0.9290 0.6940 0.1250]; p.MarkerEdgeColor = [0.9290 0.6940 0.1250]; 
p.MarkerSize = 10;

text(x(end),massEffFries,'\leftarrow Fries etl 2017')

%%  Koc etl. 2014 http://dx.doi.org/10.1016/j.trb.2014.09.008
curbWeightKoc = [3500 5500 14000]; %14 ton
payLoadKoc = [4000 12500 26000]; 
grossWeightKoc = curbWeightKoc + payLoadKoc; 
massEffKoc = 100 - ((curbWeightKoc./grossWeightKoc)*100); 

p = plot(x(end-1),massEffKoc(1),"pentagram",MarkerSize=10); 
p.MarkerFaceColor = [0.9290 0.6940 0.1250]; p.MarkerEdgeColor = [0.9290 0.6940 0.1250]; 
text(x(end-1),massEffKoc(1),'\leftarrow M. Koc etl. 2014')

p = plot(x(end-1),massEffKoc(2),"pentagram",MarkerSize=10); 
p.MarkerFaceColor = [0.9290 0.6940 0.1250]; p.MarkerEdgeColor = [0.9290 0.6940 0.1250]; 
text(x(end-1),massEffKoc(2),'\leftarrow M. Koc etl. 2014')

p = plot(x(end),massEffKoc(3),"pentagram",MarkerSize=10); 
p.MarkerFaceColor = [0.9290 0.6940 0.1250]; p.MarkerEdgeColor = [0.9290 0.6940 0.1250]; 
text(x(end),massEffKoc(3),'\leftarrow M. Koc etl. 2014')

%% Payload to Structural mass ISBN 978-0-309-49635-3 | DOI 10.17226/25542
% The weight distribution of major components catergories in class 8 truck!
% Chassis + frame 12 % 
% truck body structure 19 %
% So here im intrested in investigating how much structural weight are
% needed to carry the payload? 

structuralMassGhandariz = grossWeightGhandariz .* 0.32;
payloadEffGhandariz = (payLoadGhandariz./structuralMassGhandariz)*100

structuralMassFries = grossWeightFries .* 0.32;
payloadEffFries = (payLoadFries./ structuralMassFries) .* 100
 
structuralMassKoc = grossWeightKoc .* 0.32; 
payloadEffKoc = (payLoadKoc ./ structuralMassKoc).*100



%% EV trucks 
% volta zero Electric trucks https://voltatrucks.com/volta-zero-ambient
voltaTrucksPaylodCap = [7650 7150]; 
voltaTrucksGVW = [16000 16000]; 
% Nikola trucks http://www.nikolamotor.com/wp-content/uploads/2023/04/BEV028_TreBEV-US-6x2-Spec-Sheet-01.13.2023.pdf
NikolaTruckCurbWeight = 29297; % lb
NikolaTruckGVW = 80000; % lb


%% Illustration 


studyNames = {'Ghandariz etl. 2021' 'Ghandariz etl. 2021' 'M. Fries etl 2017' 'Koc etl. 2014' 'Koc etl. 2014' 'Koc etl. 2014','Volta Zero S R','VoltaZero LR','TRE BEV',''}; 
dataPints = [payloadEffGhandariz payloadEffFries payloadEffKoc]; 
figure2 = figure(2); clf, 
hold on 
x = linspace(1,numel(dataPints),numel(dataPints)); 
p = plot(x,dataPints,"hexagram",MarkerSize=20); 
p.MarkerFaceColor = [0.6350 0.0780 0.1840];
p.MarkerEdgeColor = [0.6350 0.0780 0.1840] 
axis([0 7 150 230]); 

ylabel('[%]')
xticks(x)
xticklabels(studyNames)
grid on
ax = gca; 
ax.TickLabelInterpreter = "latex";
ax.FontSize = 20;
pictureWidth = 21; % This one is picture width
hw_ratio = 0.8;  
set(figure1,'Units','centimeters','Position',[60 50 pictureWidth hw_ratio*pictureWidth]);
[t,s] = title('paylod-to-structural mass ratio','FontWeight','bold','Interpreter','latex');





clear all 
clc 
%% DriveCycle profile 
h = 1; % duraion of the driveCycle 
nDriveCyceles = 50; 
avrageSpeed = linspace(50,90,nDriveCyceles); 
Ttot = 60*60*h; 
timeSeries = linspace(0,Ttot,Ttot+1); 
driveCycleNames = ''; 
velFirstSectionIndex = ceil(numel(timeSeries).* 0.10); 
velMidsectionIndex   =   ceil(numel(timeSeries).* 0.9); 
velprofile = zeros(nDriveCyceles,Ttot);  
timeSeriesProfile = zeros(nDriveCyceles,Ttot+1);


%% DrivreCycle generation 

for i = 1:nDriveCyceles
velprofile(i,1:velFirstSectionIndex) = linspace(0,avrageSpeed(i),velFirstSectionIndex);
velprofile(i,velFirstSectionIndex+1:velMidsectionIndex) = avrageSpeed(i); 
zeroIndexes = find(~velprofile(i,:)); 
velprofile(i,velMidsectionIndex+1:end) = linspace(avrageSpeed(i),0,numel(zeroIndexes)-1);
timeSeriesProfile(i,:) = timeSeries; 
end 

timeSeriesProfile(:,end) = []; 

t = timeSeriesProfile; 
vel = velprofile; 


%% save the drive cycles with diferent drive cycle. 

save("driveCyclesLonghaul.mat","t","vel","-v7.3")

clear all 
clc
%% 
global model xnam 
addpath('./truckModel')
addpath('./driveCycle')
import truckModel.* driveCycle.* 
%% Initiate the truck model again
model = initruckModel;
model.objfunc = 'TransportEff';
model = updateDependentVars(model);
model.drivecycle = driveCycle.getGeneratedDriveCycleData('driveCyclesLonghaul');
globalDriveCycles = driveCycle; 
globalDriveCycles = driveCycle.setGlobalDriveCycle(globalDriveCycles,model.drivecycle); 
tDemand = [17.1449   18.9841   20.8150   22.6500   24.4810   26.3078   28.1510   29.9779   31.8294   33.6397   35.5159   37.3756]; 
%% Set optimisatoin parameters
xval=[model.beam.scale model.beam.L model.load.q(2).capacity]';
xnam={'beam.scale' 'beam.L' 'load.q(2).capacity' }';
xmin=[0.1 0.1 0.1]';
xmax=[10 50 150]';
%% pso 
% options = optimoptions('particleswarm','Display','iter','FunctionTolerance', 1e-2);
% nvars = 3; 
% fun = @getobjectiveValuePSO; 
    
%% fmincon
fun = @objval;
nlcon = @nonlcon;
options = optimoptions('fmincon','Display','iter');
%% Initiate optimisation runs 
[n,m] = size(tDemand); 
resultsOfAlteringTdemand =  cell(n,4); 
for i = 1:n
    globalDriveCycles.currentSimulation = 50;
    model.drivecycle = driveCycle.setDriveCycleforSimulation(globalDriveCycles); 
    model.mission.demand = tDemand; 
    %[x,fval,exitflag,output] = particleswarm(fun,nvars,xmin.',xmax.',options)
    [x,fval,exitflag,output] = fmincon(fun,xval,[],[],[],[],xmin.',xmax.',nlcon,options);
    resultsOfAlteringTdemand{i,1} = model;
    resultsOfAlteringTdemand{i,2} = x;
    resultsOfAlteringTdemand{i,3} = fval;
    resultsOfAlteringTdemand{i,4} = output;
	i
end


%% Save data
save("Results\alteringTdemand.mat","resultsOfAlteringTdemand","-v7.3")

%% Fmincon
function F = objval(x)
global model xnam 
for i=1:numel(x)
    eval(['model.',xnam{i},'=x(i);'])
end
model=truckModel.updateDependentVars(model);
F = truckModel.computeEnergy(model).*1e-7;

end

function [Cineq,Ceq] = nonlcon(x)
global model xnam

for i=1:numel(x)
    eval(['model.',xnam{i},'=x(i);'])
    model = truckModel.updateDependentVars(model);
    deflection = max(abs(model.beam.delta(linspace(0,model.beam.L,1e4))));
    Cineq = [(deflection - model.fmax) (model.mission.demand-model.load.q(2).capacity)] ;
    Ceq = [];
end


end


%% Particle swarn 
function F = getobjectiveValuePSO(x)
global model xnam 
for i=1:numel(x)
    eval(['model.',xnam{i},'=x(i);'])
end
model = truckModel.updateDependentVars(model);
totalEnergyConsumption = truckModel.computeEnergy(model).*1e-7;
%energyPerkm = totalEnergyConsumption / (model.drivecycle.totalx / 1e3) ;  

F = (model.mission.Ndelivery .* totalEnergyConsumption)   + (model.penaltyVal.* (getMaximumDeflectionValue(model)> model.fmax)); 

end

function fval = getMaximumDeflectionValue(model)
deflection = max(abs(model.beam.delta(linspace(0,model.beam.L,1e4))));
fval(1) = deflection;
end
% This file is part of ecoOptimize, a code to optimize a design model for
% minimum eco impacts subject to functional requirements.
%
% Copyright (C) 2020 Ciarán O'Reilly <ciaran@kth.se> &
% Khashayar Shahrezaei <khasha@kth.se>
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.%
%
clear all

%%
global model history xnam 
addpath('./truckModel')
addpath('./driveCycle')
import truckModel.* driveCycle.* dataFit.*
%% Initiate the truck model again
solver = 'PSO'; % GCMMA PSO 
model = initruckModel;
model.objfunc = 'TransportEff';
model = updateDependentVars(model);
model.drivecycle = driveCycle.cycle('REGIONAL3');

%% Set optimisatoin parameters PSO
% xval=[model.beam.scale model.beam.L model.load.q(2).capacity model.mission.Ndelivery]';
% xnam={'beam.scale' 'beam.L' 'load.q(2).capacity' 'mission.Ndelivery'}';
% xmin=[0.1 0.1 0.1 0.1]';
% xmax=[10 50 600 10]';

%% Set optimisatoin parameters PSO
xval=[model.beam.scale model.beam.L model.load.q(2).capacity ]';
xnam={'beam.scale' 'beam.L' 'load.q(2).capacity'}';
xmin=[0.1 1 50]';
xmax=[10 30 150]';

%% Initiate optimisation
switch solver
    case 'GCMMA'
        gcmma=GCMMA.init(@optFuncs,xval,xnam,xmin,xmax);
        disp(['Optimizing for: ',model.objfunc])
        gcmma.displive=1;
        figure(2), clf, gcmma.plotlive=1;
        gcmma.maxoutit=20;
        [gcmma,xval]=GCMMA.run(gcmma);
        [f0val,fval]=optFuncs(xval,xnam,false);
    case 'PSO'
        figure(4); 
        clf; 
        options = optimoptions('particleswarm','Display','iter','FunctionTolerance', 1e-2,'OutputFcn',@plotValues);
        nvars = 3; 
        fun = @getobjectiveValuePSO; 
        [x,fval,exitflag,output] = particleswarm(fun,nvars,xmin.',xmax.',options)


    otherwise
        history = [];
        fun = @objval;
        nlcon = @nonlcon;
        options = optimoptions('fmincon','Display','iter');

        [x,fval,exitflag,output] = fmincon(fun,xval,[],[],[],[],xmin.',xmax.',nlcon,options);
end
%% Post process
figure(3), clf, truckModel.dispBeamFuncs(model.beam)
X2 = (model.beam.L+0.5);
Y2 = (model.shape.H +0.5);
figure(2),clf,truckModel.dispTruck(model)
axis([-0.5 X2,-0.5 Y2])

fuelEcoLoaded   = sum(truckModel.computeEnergy(model,0)) * ceil(model.mission.Ndelivery);
fuelEcoEmpty    = sum(truckModel.computeEnergy(model,1)) * ceil(model.mission.Ndelivery);
fuelEco = fuelEcoLoaded + fuelEcoEmpty;

vehTonneKilometer   = model.load.q(2).capacity * ceil(model.mission.Ndelivery)*2* (model.drivecycle.cr/(1000)/ceil(model.mission.Ndelivery)*2)
energyPerPayload    = fuelEco/ (model.load.q(2).capacity *  ceil(model.mission.Ndelivery))

%% Fmincon
function F = objval(x)
global model xnam 
for i=1:numel(x)
    eval(['model.',xnam{i},'=x(i);'])
end
model=truckModel.updateDependentVars(model);
totalEnergyConsumption = truckModel.computeEnergy(model).*1e-7;
energyPerkm = totalEnergyConsumption / (model.drivecycle.totalx / 1e3) ;  

F = totalEnergyConsumption;  

end


function [Cineq,Ceq] = nonlcon(x)
global model xnam

for i=1:numel(x)
    eval(['model.',xnam{i},'=x(i);'])
    model = truckModel.updateDependentVars(model);
    deflection = max(abs(model.beam.delta(linspace(0,model.beam.L,1e4))));
    Cineq = [deflection - model.fmax] ;
    Ceq = [];
end


end

function stop = myoutput(x,optimvalues,state)
global history effMeasures
stop = false;
if isequal(state,'iter')

   history = [history; x' optimvalues.fval effMeasures.mass effMeasures.shape effMeasures.service effMeasures.productivity];

end
end

%% Particle swarn 
function F = getobjectiveValuePSO(x)
global model xnam 
for i=1:numel(x)
    eval(['model.',xnam{i},'=x(i);'])
end
model = truckModel.updateDependentVars(model);
totalEnergyConsumption = truckModel.computeEnergy(model).*1e-7 .* model.mission.Ndelivery;
%TotalCapacity = model.load.q(2).capacity * model.mission.Ndelivery; 
%energyPerpayloadPerkm = (totalEnergyConsumption /TotalCapacity) / (model.drivecycle.totalx / 1e3) ;  
F = (totalEnergyConsumption) + (model.penaltyVal.* (getMaximumDeflectionValue(model) > model.fmax))  + (model.penaltyVal.* ((ceil(model.mission.Ndelivery) .* model.load.q(2).capacity) < model.mission.demand)) ; 
end

function fval = getMaximumDeflectionValue(model)
deflection = max(abs(model.beam.delta(linspace(0,model.beam.L,1e4))));
fval(1) = deflection;
end 

function stop = plotValues(optimValues,state)
global xnam
stop = false; % This function does not stop the solver
switch state
    case 'init'
        nplot = size(optimValues.swarm,2); % Number of dimensions
        for i = 1:nplot % Set up axes for plot
            subplot(nplot,1,i);
            tag = sprintf('psoplotrange_var_%g',i); % Set a tag for the subplot
            plot(optimValues.iteration,0,'-*k','Tag',tag); 
            ylabel(char(xnam(i,1)));
        end
        xlabel('Iteration','interp','none'); % Iteration number at the bottom
        subplot(nplot,1,1) % Title at the top
        title('Log range of particles by component')
        setappdata(gcf,'t0',tic); % Set up a timer to plot only when needed
    case 'iter'
        nplot = size(optimValues.swarm,2); % Number of dimensions
        for i = 1:nplot
            subplot(nplot,1,i);
            % Calculate the range of the particles at dimension i
            irange = max(optimValues.swarm(:,i)) - min(optimValues.swarm(:,i));
            tag = sprintf('psoplotrange_var_%g',i);
            plotHandle = findobj(get(gca,'Children'),'Tag',tag); % Get the subplot
            xdata = plotHandle.XData; % Get the X data from the plot
            newX = [xdata optimValues.iteration]; % Add the new iteration
            plotHandle.XData = newX; % Put the X data into the plot
            ydata = plotHandle.YData; % Get the Y data from the plot
            newY = [ydata irange]; % Add the new value
            plotHandle.YData = newY; % Put the Y data into the plot
        end
        if toc(getappdata(gcf,'t0')) > 1/30 % If 1/30 s has passed
            drawnow % Show the plot
            setappdata(gcf,'t0',tic); % Reset the timer
        end
    case 'done'
        % No cleanup necessary
end
end







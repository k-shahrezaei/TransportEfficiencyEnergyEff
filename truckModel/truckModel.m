classdef truckModel
    
    properties
        obj
    end
    
    methods(Static)
                
        function obj = initruckModel
            obj.beam.L = 10;
            obj.beam.scale = 1.0;
            obj.beam.H0 = [0.04 0.4 0.04];
            obj.beam.B0 = [0.4 0.08 0.4];
            
            obj.beam.E = 2.1e11;
            obj.beam.rho = 7800;
            
            obj.axle.rad = 0.4;
            obj.axle.Crr = 0.01;
            obj.axle.FRV1= 4978200;      % [kN/m]  Tyre stiffness, vertival tire rate
            obj.axle.m1  =0.44;          % [m] suspension roll axis height
            obj.axle.CDG1=184000;        % [kNm/rad] suspension roll stiffness for roll axis
            obj.axle.mTire = 2;        % [kg]
            
            obj.shape.B = 2.0;      % [m] nominal truck width
            obj.shape.Cd = 0.6;
            obj.shape.CdFunction = truckModel.CdSurrogateFunction; 

            
            obj.load.q(2).capacity = 100; % initial loading capacity
            obj.load.q(2).rho = 500;
            obj.load.q(2).mass = obj.load.q(2).rho .* obj.load.q(2).capacity;

            obj.load.q(3).mass  = 3000; % mass of the battery  
            obj.load.q(3).rho = 1000;   % density of battery pack

            obj.load.q(4).mass = 2 * 165;
            obj.load.q(4).rho = 1000;       % density of elerctric motor
            
            obj.load.q(5).mass = 2 * 165;  
            obj.load.q(5).rho = 1000;       % density of elerctric motor
            
            
            obj.mission.demand = 113; % the daily demand [m^3]
            obj.mission.Ndelivery = []; 

            
            
            obj.fmax= 1e-2;
            obj.penaltyVal = 1e3; 
            

        end
        
        function obj=updateDependentVars(obj)
            obj.beam.B = obj.beam.B0*obj.beam.scale;
            obj.beam.H = obj.beam.H0*obj.beam.scale;
            I0=obj.beam.B.*obj.beam.H.^3/12;
            d=([0 cumsum(obj.beam.H(1:end-1))]+obj.beam.H/2-sum(obj.beam.H)/2);
            A=obj.beam.B.*obj.beam.H;
            I=I0+A.*d.^2;
            obj.beam.A=sum(A);
            obj.beam.I=sum(I);
            obj.load.q(1).x = [0 obj.beam.L];
            obj.load.q(1).fun = str2func(['@(x)',mat2str(-obj.beam.rho*obj.beam.A*9.81),'+0.*x']);
            % here you need to update the new loading mass%%%%%%%%%%%%%%%%%%%%%
            obj.load.q(2).mass = obj.load.q(2).rho .* obj.load.q(2).capacity;%%
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            obj.mission.Ndelivery = ceil(obj.mission.demand ./ ( obj.load.q(2).capacity));
            obj.load.q(2).x = [0 obj.beam.L];
            obj.load.q(3).x = [0.1*obj.beam.L 0.5*obj.beam.L]; 
            obj.load.q(4).x = [0.05*obj.beam.L   0.15*obj.beam.L]; 
            obj.load.q(5).x = [0.75*obj.beam.L   0.85*obj.beam.L]; 
            
            
            for i = 2:numel(obj.load.q)
                obj.load.q(i).fun = str2func(['@(x)',mat2str(-obj.load.q(i).mass/abs(diff(obj.load.q(i).x))*9.81),'+0.*x']);
                obj.load.q(i).H = obj.load.q(i).mass/obj.load.q(i).rho/abs(diff(obj.load.q(i).x))/obj.shape.B;
                if i == 2
                    obj.load.q(i).x(2) = obj.load.q(i).mass/obj.load.q(i).rho/obj.load.q(i).H/obj.shape.B;
                end
            end
            obj.load.H = max([obj.load.q.H]);
            obj.axle.R(1).x = 0.1*obj.beam.L;
            obj.axle.R(2).x = 0.8*obj.beam.L;
            obj.shape.H = 2 * obj.axle.rad + sum(obj.beam.H) + obj.load.H;
            obj.beam.L = obj.load.q(2).mass/obj.load.q(2).rho/obj.load.q(2).H/obj.shape.B;
            obj.shape.L = obj.beam.L;
            obj.shape.A = obj.shape.H*obj.shape.B;
            obj = truckModel.computeReactions(obj);
            obj = truckModel.evaluateBeamFuncs(obj);             
        end
        
        function obj = computeReactions(obj)
            F = arrayfun(@(y) integral(y.fun,y.x(1),y.x(2)),obj.load.q);
            M = arrayfun(@(y) integral(@(x)y.fun(x).*x,y.x(1),y.x(2)),obj.load.q);
            r = [obj.axle.R(:).x];
            invA = [-r(2)/(r(1) - r(2)),  1/(r(1) - r(2))
                r(1)/(r(1) - r(2)), -1/(r(1) - r(2))];
            a = -invA*[sum(F);sum(M)];
            obj.axle.R(1).a = a(1);
            obj.axle.R(2).a = a(2);
        end
        
        function obj = evaluateBeamFuncs(obj)
            Q=sum([arrayfun(@(y) y.a*heaviside(sym('x')-y.x),obj.axle.R) arrayfun(@(y) int(sym(y.fun),y.x(1),sym('x'))*heaviside(sym('x')-y.x(1)),obj.load.q) arrayfun(@(y) -int(sym(y.fun),y.x(2),sym('x'))*heaviside(sym('x')-y.x(2)),obj.load.q)]);
            M=int(Q,0,sym('x'));
            T=(int(M,0,sym('x'))-int(int(M,0,sym('x')),obj.axle.R(1).x,obj.axle.R(2).x)/(obj.axle.R(2).x-obj.axle.R(1).x))/obj.beam.E/obj.beam.I;
            W=int(T,obj.axle.R(1).x,sym('x'));
            obj.beam.shear = matlabFunction(Q);
            obj.beam.moment = matlabFunction(M);
            obj.beam.curve = matlabFunction(T);
            obj.beam.delta = matlabFunction(W);
        end
        
        
        
        function mass=computeMass(model)
            model=model.beam;
            mass=model.A.*model.L.*model.rho/1e3; %scale so it's not too large
        end
        
        
        function Energy = computeEnergy(model,empthy)
            if nargin < 2
                empthy = 0;
            end 
            crr=model.axle.Crr*(1-model.drivecycle.r); %[-] %coefficient of effective rolling resistance
            cd = model.shape.CdFunction(model.beam.L / model.shape.H);  %[-] %coefficient of drag
            CR=model.drivecycle.cr; %[m] drive cycle rolling resistance constant
            CA=model.drivecycle.ca; %[m^2/s^2] drive cycle acceleration constant
            CD=model.drivecycle.cd; %[m^3/s^2] drive cycle aerodynamic constant
            diffEff=0.42; %[-] differential efficiency (petrol)
            if empthy
                mass = truckModel.computeMass(model)*1e3; % Total mass (Beam + Payload)
            else
                mass = truckModel.computeMass(model)*1e3 + sum([model.load.q(2:end).mass]); % Total mass (Beam + Payload)
            end
            A=model.shape.A; %total Frontal area
            
            EA_km=(CA)/(CR/1e3)/diffEff.* mass; %[J/100km]
            ER_km=(9.81*crr*CR)/(CR/1e3)/diffEff.* mass; %[J/100km]
            ED_km=(0.5*1.2*cd*CD)/(CR/1e3)/diffEff.*A; %[J/100km]
            Energy= sum([EA_km ER_km ED_km]); %[J/100km]

            % EA=(CA)/diffEff.* mass; %[J]
            % ER=(9.81*crr*CR)/diffEff.* mass; %[J]
            % ED=(0.5*1.2*cd*CD)/diffEff.*A; %[J]
            % Energy= sum([EA ER ED]); %[J]
        end
        function dispBeamFuncs(beam)
            subplot(2,2,1)
            fplot(beam.shear,[0 beam.L])
            grid on, xlabel('x'),ylabel('shear force')
            subplot(2,2,2)
            fplot(beam.moment,[0 beam.L])
            grid on, xlabel('x'),ylabel('bending moment')
            subplot(2,2,3)
            fplot(beam.curve,[0 beam.L])
            grid on, xlabel('x'),ylabel('curvature')
            subplot(2,2,4)
            fplot(beam.delta,[0 beam.L])
            grid on, xlabel('x'),ylabel('delta')
        end
        function dispBeamSection(model,fill)
            N=max([numel(model.B) numel(model.H)]);
            B=repmat(model.B,1,N-numel(model.B)+1);
            H=repmat(model.H,1,N-numel(model.H)+1);
            H0=([0 cumsum(H(1:end-1))]-sum(H)/2);
            B0=-B/2;
            C=colormap;
            for i=1:N
                R=rectangle('Position',[B0(i) H0(i) B(i) H(i)]);
                if fill
                    set(R,'Facecolor',C(50*(i-1)+1,:))
                end
            end
            axis([B0(i)-B(i)/10 B0(i)+B(i)+B(i)/10 H0(i)-H(i)/10 H0(i)+H(i)+H(i)/10])
            axis auto, axis equal
            xlabel('b [m]'), ylabel('h [m]')
        end
        
        function dispTruck(obj)
            C=colormap;
            r=obj.axle.rad;
            for i = 1:2
                x = obj.axle.R(i).x;
                R=rectangle('Position',[x-r 0 2*r 2*r],'Curvature',[1 1]);
                set(R,'Facecolor','k')
            end
            h0 = 2*r;
            for i = 1:3
                h = obj.beam.H(i);
                L = obj.beam.L;
                R=rectangle('Position',[0 h0 L h]);
                set(R,'Facecolor','g')
                h0 = h0 + h;
            end
            for i = 2:numel(obj.load.q)
                R=rectangle('Position',[obj.load.q(i).x(1) h0 obj.load.q(i).x(2)-obj.load.q(i).x(1) obj.load.q(i).H]);
                set(R,'Facecolor',C(50*(i-1)+1,:))
            end
            rectangle('Position',[0 r obj.beam.L obj.shape.H],'LineStyle','--')
            axis equal
            xlabel('l [m]'), ylabel('h [m]')
        end


        function fitresult = CdSurrogateFunction 
            aspectsRatio =  [0.5 1 2 3 4 6 10];
            CD =            [1 0.9 0.8 0.82 0.85 0.9 1];


            [xData, yData] = prepareCurveData( aspectsRatio, CD );

            % Set up fittype and options.
            ft = fittype( 'exp2' );
            opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
            opts.Display = 'Off';
            opts.StartPoint = [1.40258995651883 -3.52080425672074 0.747461440667134 0.0301297152039899];

            % Fit model to data.
            [fitresult, gof] = fit( xData, yData, ft, opts );
            

        end
        
        
    end
end
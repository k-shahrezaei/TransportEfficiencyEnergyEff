classdef driveCycle

    properties
        currentSimulation
        totalx
        scaleTime
        t
        x
        vel
        avgVel 
        acc
        r
        cr
        cd
        ca
    end

    methods(Static)

        function obj = cycle(drivecycle)
            if ~isstruct(drivecycle)
                name = drivecycle;
                drivecycle = [];
                drivecycle.name = name;
            end
            tabledata = (drivecycle.name + ".xlsx");
            obj.name = name;
            t = readtable(char(tabledata));
            vel = t.vel/3.6; %km/h to m/s
            t = t.t;
            n=length(t);
            D=spdiags(repmat([-0.5 0 0.5],n,1),-1:1,n,n)+sparse([1 1 n n],[1 2 n-1 n],[-1 0.5 -0.5 1]);
            dt=D*t;
            dx=vel.*dt;
            x=cumsum(dx);
            acc=D*vel./dt;
            r=sum(acc<0)/numel(acc);
            acc(acc<0)=0;
            obj.totalx = sum(x);
            obj.scaleTime = 1;
            obj.t=t;
            obj.x=x;
            obj.vel=vel;
            obj.avgVel = mean(vel);
            obj.acc=acc;
            obj.r=r;
            obj.cr=sum(dx);
            obj.cd=sum(vel.^3.*dt);
            obj.ca=sum(acc.*vel.*dt);
        end

        function obj = getGeneratedDriveCycleData(drirveCycleProfile)
            if ~isstruct(drirveCycleProfile)
                name = drirveCycleProfile;
                drivecycle = [];
                drivecycle.name = name;
            end
            tabledata = (drivecycle.name + ".");
            obj.name = name;
            t = load(char(tabledata));
            vel = t.vel/3.6; %km/h to m/s
            t = t.t;
            [m,n]=size(t);
            D=spdiags(repmat([-0.5 0 0.5],n,1),-1:1,n,n)+sparse([1 1 n n],[1 2 n-1 n],[-1 0.5 -0.5 1]);
            dt=D*t(1,:).';
            dx = zeros(m,n);
            x = zeros(m,n);
            acc = zeros(m,n);
            r = zeros(m,1);
            for i = 1:m
                dx(i,:)=vel(i,:).*dt.';
                x(i,:)=cumsum(dx(i,:));
                acc(i,:)=(D*vel(i,:).') ./ dt;
                r(i,:)=sum(acc(i,:)<0)/numel(acc(i,:));
            end
            acc(acc<0)=0;
            obj.totalx = sum(x,2);
            obj.scaleTime = 1;
            obj.t=t;
            obj.x=x;
            obj.vel=vel;
            obj.avgVel = mean(vel,2);
            obj.acc=acc;
            obj.r=r;
            obj.cr=sum(dx,2);
            obj.cd=sum(vel.^3.* dt.',2);
            obj.ca=sum(acc.*vel.*dt.',2);
            
        end

        function obj = setDriveCycleforSimulation(currentDriveCycle)
            k   = currentDriveCycle.currentSimulation;
            obj.totalx =  currentDriveCycle.totalx(k); 
            obj.scaleTime = currentDriveCycle.scaleTime;
            obj.t=currentDriveCycle.t(k,:); 
            obj.x=currentDriveCycle.x(k,:); 
            obj.vel=currentDriveCycle.vel(k,:); 
            obj.avgVel = currentDriveCycle.avgVel(k); 
            obj.acc =currentDriveCycle.acc(k,:); 
            obj.r=currentDriveCycle.r(k,:); 
            obj.cr=currentDriveCycle.cr(k,:); 
            obj.cd=currentDriveCycle.cd(k,:); 
            obj.ca=currentDriveCycle.ca(k,:); 
        end
        function globalDriveCycles = setGlobalDriveCycle(globalDriveCycles,driveCycles)
            globalDriveCycles.totalx =  driveCycles.totalx; 
            globalDriveCycles.scaleTime = driveCycles.scaleTime;
            globalDriveCycles.t=driveCycles.t; 
            globalDriveCycles.x=driveCycles.x; 
            globalDriveCycles.vel=driveCycles.vel; 
            globalDriveCycles.avgVel = driveCycles.avgVel; 
            globalDriveCycles.acc =driveCycles.acc; 
            globalDriveCycles.r=driveCycles.r; 
            globalDriveCycles.cr=driveCycles.cr; 
            globalDriveCycles.cd=driveCycles.cd; 
            globalDriveCycles.ca=driveCycles.ca; 
        
        end 


    end
end


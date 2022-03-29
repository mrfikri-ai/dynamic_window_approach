% DWA Pseudocode
% BEGIN DWA(robotPose,robotGoal,robotModel)  
%    laserscan = readScanner()  
%    allowable_v = generateWindow(robotV, robotModel)  
%    allowable_w  = generateWindow(robotW, robotModel)  
%    for each v in allowable_v  
%       for each w in allowable_w  
%       dist = find_dist(v,w,laserscan,robotModel)  
%       breakDist = calculateBreakingDistance(v)  
%       if (dist > breakDist)  //can stop in time  
%          heading = hDiff(robotPose,goalPose, v,w)   
%           //clearance
%          clearance = (dist-breakDist)/(dmax - breakDist)   
%          cost = costFunction(heading,clearance, abs(desired_v - v))  
%          if (cost > optimal)  
%             best_v = v  
%             best_w = w  
%             optimal = cost  
%     set robot trajectory to best_v, best_w  
% END


function [] = DynamicWindowApproachSample()
 
close all;
clear all;
 
disp('Dynamic Window Approach sample program start!!')

% Initial state of the robot [x (m), y (m), yaw (Rad), v (m / s), w (rad / s)]
x=[0 0 pi/2 0 0]';

% Target point position [x (m), y (m)]
goal=[7,8];

% Obstacle position list [x (m) y (m)]
obstacle=[0 2;
          4 2;
          4 4;
          5 4;
          5 5;
          5 6;
          5 9
          8 8
          8 9
          7 9
          6 5
          6 3
          6 8
          6 7
          7 4
          9 8
          9 11
          9 6];
      
% Obstacle radius for collision detection
obstacleR=0.3;

% Time [s]
global dt; dt=0.1;

% Kinematics model
% Maximum speed m / s], maximum rotation speed [rad / s], acceleration [m / ss], rotation acceleration [rad / ss],
% Speed resolution rate [m / s], rotation speed resolution rate [rad / s]]
Kinematic=[1.0,toRadian(20.0),0.2,toRadian(50.0),0.01,toRadian(1)];

% Parameter number reference [heading, dist, velocity, predictDT]
evalParam=[0.05,0.2,0.1,3.0];

% area of the environment
area=[-1 11 -1 11];

% Imitation experiment
result.x=[];
tic;
% movcount=0;
% Main loop
for i=1:5000
    % DWA 
    [u,traj]=DynamicWindowApproach(x,Kinematic,goal,evalParam,obstacle,obstacleR);
    x=f(x,u); % The robot moves to the next moment
    
    % save the simulation results
    result.x=[result.x; x'];
    
    % when reach the destination
    if norm(x(1:2)-goal')<0.5
        disp('Arrive Goal!!');break;
    end
    
    %====Animation====
    hold off;
    ArrowLength=0.5;% 
    
    % Robot
    quiver(x(1),x(2),ArrowLength*cos(x(3)),ArrowLength*sin(x(3)),'ok');hold on;
    plot(result.x(:,1),result.x(:,2),'-b');hold on;
    plot(goal(1),goal(2),'*r');hold on;
    plot(obstacle(:,1),obstacle(:,2),'*k');hold on;
    
    % Explore the track
    if ~isempty(traj)
        for it=1:length(traj(:,1))/5
            ind=1+(it-1)*5;
            plot(traj(ind,:),traj(ind+1,:),'-g');hold on;
        end
    end
    axis(area);
    grid on;
    drawnow;
    %movcount=movcount+1;
    %mov(movcount) = getframe(gcf);% 
end
toc
%movie2avi(mov,'movie.avi');
 

function [u,trajDB]=DynamicWindowApproach(x,model,goal,evalParam,ob,R)% DWA²ÎÊýÊäÈë

% Dynamic Window [vmin,vmax,wmin,wmax]
Vr=CalcDynamicWindow(x,model);

% Numerical calculation
[evalDB,trajDB]=Evaluation(x,Vr,goal,ob,R,model,evalParam);

if isempty(evalDB)
    disp('no path to goal!!');
    u=[0;0];return;
end

% Regularization of each function
evalDB=NormalizeEval(evalDB);

% Calculation of the final evaluation function
feval=[];
for id=1:length(evalDB(:,1))
    feval=[feval;evalParam(1:3)*evalDB(id,3:5)'];
end
evalDB=[evalDB feval];

[maxv,ind]=max(feval);% optimal evaluation function
u=evalDB(ind,1:2)';% 

function [evalDB,trajDB]=Evaluation(x,Vr,goal,ob,R,model,evalParam)
%
evalDB=[];
trajDB=[];
for vt=Vr(1):model(5):Vr(2)
    for ot=Vr(3):model(6):Vr(4)
        % Trajectory estimation; get xt: the predicted pose after 
        % the robot moves forward; traj: the trajectory between the current moment and the predicted moment
        
        %evalParam(4), forward simulation time;
        [xt,traj]=GenerateTrajectory(x,vt,ot,evalParam(4),model);  
        
        % Calculation of evaluation
        heading=CalcHeadingEval(xt,goal);
        dist=CalcDistEval(xt,ob,R);
        vel=abs(vt);
        
        % braking system
        stopDist=CalcBreakingDist(vel,model);
        if dist>stopDist % 
            evalDB=[evalDB;[vt ot heading dist vel]];
            trajDB=[trajDB;traj];
        end
    end
end

function EvalDB=NormalizeEval(EvalDB)
% Regularization of the function
if sum(EvalDB(:,3))~=0
    EvalDB(:,3)=EvalDB(:,3)/sum(EvalDB(:,3));
end
if sum(EvalDB(:,4))~=0
    EvalDB(:,4)=EvalDB(:,4)/sum(EvalDB(:,4));
end
if sum(EvalDB(:,5))~=0
    EvalDB(:,5)=EvalDB(:,5)/sum(EvalDB(:,5));
end

function [x,traj]=GenerateTrajectory(x,vt,ot,evaldt,model)
% Trajectory generation function
% evaldt: forward simulation time; vt, ot current velocity and angular velocity;

global dt;
time=0;
u=[vt;ot]; % input value
traj=x; % robot trajectory
while time<=evaldt
    time=time+dt; % time update
    x=f(x,u); % motion update
    traj=[traj x];
end

function stopDist=CalcBreakingDist(vel,model)
% Calculate the braking distance according to kinematic model
global dt;
stopDist=0;
while vel>0
    stopDist=stopDist+vel*dt; % Stop or braking distance calculation
    vel=vel-model(3)*dt;% 
end

function dist=CalcDistEval(x,ob,R)
% Obstacle distance evaluation

dist=100;
for io=1:length(ob(:,1))
    disttmp=norm(ob(io,:)-x(1:2)')-R;
    if dist>disttmp % Minimum distance to obstacle
        dist=disttmp;
    end
end

% The obstacle distance evaluation limits a maximum value. 
%If it is not set, once a trajectory has no obstacles, it will take too much weight
if dist>=2*R
    dist=2*R;
end

function heading=CalcHeadingEval(x,goal)
% evaluation of heading robot

theta=toDegree(x(3)); % Robot orientation

% Orientation of target point
goalTheta=toDegree(atan2(goal(2)-x(2),goal(1)-x(1))); 

if goalTheta>theta
    targetTheta=goalTheta-theta;% [deg]
else
    targetTheta=theta-goalTheta;% [deg]
end

heading=180-targetTheta;

function Vr=CalcDynamicWindow(x,model)
%
global dt;
% Maximum and minimum range of car speed
Vs=[0 model(1) -model(2) model(2)];

% Dynamic window calculated based on current speed and acceleration limits
Vd=[x(4)-model(3)*dt x(4)+model(3)*dt x(5)-model(4)*dt x(5)+model(4)*dt];

% Final Dynamic Window
Vtmp=[Vs;Vd];
Vr=[max(Vtmp(:,1)) min(Vtmp(:,2)) max(Vtmp(:,3)) min(Vtmp(:,4))];

function x = f(x, u)
% Motion Model
% u = [vt; wt];Velocity and angular velocity at the current moment
global dt;
 
F = [1 0 0 0 0
     0 1 0 0 0
     0 0 1 0 0
     0 0 0 0 0
     0 0 0 0 0];
 
B = [dt*cos(x(3)) 0
    dt*sin(x(3)) 0
    0 dt
    1 0
    0 1];

x= F*x+B*u;

function radian = toRadian(degree)
% degree to radian
radian = degree/180*pi;

function degree = toDegree(radian)
% radian to degree
degree = radian/pi*180;
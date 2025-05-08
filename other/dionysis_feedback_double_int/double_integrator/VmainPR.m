clear all; close all; clc;
% clearvars -except x y v w
clearvars -except x y v w t2 acc accy xPR yPR vPR wPR t1 accPR accyPR xN yN wN vN accN
global  a vmax vstar sigma n obs_n  q  L lambda mu2 mu1 epsilon A c p phi rho visc wmax DYNAMIC_STEP

DYNAMIC_STEP = 0;
n = 1;
obs_n = 8;
h = 0.01; A = 1; a = 5.0;
vmax = 35; wmax = 4.0; vstar = 30;

L = 6.5;
lambda = 25; % sensing radius
sigma = 5;

p = 15.0; % eccentricity of ellipse
rho = 1; % extra terms for controller, if 1 extra terms are zero
c = 2.5; % defines an area at the center of road where potential is zero
epsilon = 0.5;

t_final = 25;
 
Atol = 1*10^(-2); % the smaller the closer to continuous solution
Rtol = 1*10^(-2); % the smaller the closer to continuous solution

%% NEWTONIAN
model = 1;
visc  = 0.0;
q     = 10^(-2); % weight for potential
mu1   = 0.1;  % weight for Vstar convergence (was 0.05)
mu2   = 2.0;  % weight for lateral speed convergence

%% RELATIVISTIC
% model =0; 
% visc = 0;% 0.3*1/vmax^2;
% q = 10^(-3)/vmax^1;
% mu1 = 0.5; mu2 =1/vmax^2;
 
%  X0=x; Y0=y; V0=v; T0=w;

%% initial conditions
X0 =[5.0];
V0 =[25.0];
T0 =[-0.5];
Y0 =[3.5 - 5.0];

OBSX0 =[40.0, 40.0, 80.0, 85.0, 125.0, 120.0, 165.0, 160.0, 205.0, 200.0];
OBSV0 =[25.0, 25.0, 25.0, 25.0, 25.0, 25.0, 25.0, 25.0];
OBST0 =[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
% OBSY0 =[1.5 - 5.0, 4.0 - 5.0, 1.5- 5.0, 7.5- 5.0];
OBSY0 =[1.5 - 5.0, 7.5 - 5.0, 1.5- 5.0, 4.5- 5.0, 6.0 - 5.0, 9.0 - 5.0, 4.0 - 5.0, 1.0 - 5.0, 6.0 - 5.0, 9.0 - 5.0];

%% run
tic
% [t,xPR,yPR,wPR,vPR,accPR,accyPR,H] = heuler(h,t_final,X0,Y0,T0,V0,Atol,Rtol);
[t, xPR, yPR, wPR, vPR, obsxPR, obsyPR, obswPR, obsvPR, accPR, accyPR, Ht] = Vheuler(h, t_final, X0, Y0, T0, V0, OBSX0, OBSY0, OBST0, OBSV0, Atol, Rtol, model);
% [t,xPR,yPR,wPR,vPR,accPR,accyPR] = heun(h,t_final,X0,Y0,T0,V0);
% [t,xPR,yPR,wPR,vPR,accPR,accyPR] = Rk1(h,t_final,X0,Y0,T0,V0);
toc

if model == 1
    xN = xPR; yN = yPR; wN = wPR; vN = vPR; t1 = t; accN = accPR;
else
    xPR = xPR; yPR = yPR; wPR = wPR; vPR = vPR; t2 = t;
end

disp(max(accPR)); disp(max(accyPR)); disp(min(accPR)); disp(min(accyPR));


plot(yPR)


% max_x = max(max(xPR), max(obsxPR(:)));
% hPlot = plot(NaN, NaN, 'db');
% pbaspect([10 1 1])
% points = 100;
% mod_value = round((50 / h) / points) + 1;
% for k = 1:length(t) %loop
%     if mod(k, mod_value) == 0
%         r = rectangle('Position',[xN(k)-2.5 yN(k)-1 5 2]);
%         for o = 1:obs_n 
%             obsr(o) = rectangle('Position',[obsxPR(o,k)-2.5 obsyPR(o,k)-1 5 2]);
%         end
%         xlim([xN(k)-50,xN(k)+50])
%         ylim([-5,5])
%         pause(0.1)
%         if k < length(t) % Don't delete the very last one.
% 		    delete(r); delete(obsr);
%             % for o = 1:obs_n 
%             % end
%         end
%     end
% end

% plot(accyPR)
% plot(accPR)

%% get N values from trajectories
% points = 100;
% testAccX(points) = 0.0; testAccY(points) = 0.0;
% mod_value = round((t(end) / h) / points) + 1;
% 
% adjust = fix(length(t) / points) + 1;
% 
% testAccX(1) = accPR(1); testAccY(1) = accyPR(1); indx = 2;
% for k = 2:length(t) %loop
%     if mod(k, adjust) == 0
%         testAccX(indx) = accPR(k);
%         testAccY(indx) = accyPR(k);
%         indx = indx + 1;
%     end
% end
% disp(t(end))
% step = t(end)/(100+1);
% testX(points+1) = 0.0; testY(points+1) = 0.0; testVX(points+1) = 0.0; testVY(points+1) = 0.0;
% testX(1) = xPR(1); testY(1) = yPR(1); testVX(1) = vPR(1); testVY(1) = wPR(1);
% for k = 1 : length(testAccX)
%     testX(k + 1) = testX(k) + testVX(k)*step + 0.5*testAccX(k)*step^2;
%     testY(k + 1) = testY(k) + testVY(k)*step + 0.5*testAccY(k)*step^2;
%     testVX(k + 1) = testVX(k) + testAccX(k)*step;
%     testVY(k + 1) = testVY(k) + testAccY(k)*step;
% end
% testY = testY + 5.0;
% 
% 
% testObsX(obs_n,points+1) = 0.0; testObsY(obs_n,points+1) = 0.0; testObsVX(obs_n,points+1) = 0.0; testObsVY(obs_n,points+1) = 0.0;
% for o = 1 : obs_n
%     testObsX(o,1) = OBSX0(o); testObsY(o,1) = OBSY0(o) + 5.0; testObsVX(o,1) = OBSV0(o); testObsVY(o,1) = OBST0(o);
%     for k = 1 : length(testAccX)
%         testObsX(o,k + 1) = testObsX(o,k) + testObsVX(o,k)*0.25;
%         testObsY(o,k + 1) = testObsY(o,k) + testObsVY(o,k)*0.25;
%         testObsVX(o,k + 1) = testObsVX(o,k);
%         testObsVY(o,k + 1) = testObsVY(o,k);
%     end
% end
% 
% for k = 1:length(testY) %loop
%     % if mod(k, mod_value) == 0
%         r = rectangle('Position',[testX(k)-2.5 testY(k)-1 5 2]);
%         for o = 1:obs_n 
%             obsr(o) = rectangle('Position',[testObsX(o,k)-2.5 testObsY(o,k)-1 5 2]);
%         end
%         xlim([testX(k)-50,testX(k)+50])
%         ylim([0,10])
%         pause(0.1)
%         if k < length(t) % Don't delete the very last one.
% 		    delete(r); delete(obsr);
%             % for o = 1:obs_n 
%             % end
%         end
%     % end
% end


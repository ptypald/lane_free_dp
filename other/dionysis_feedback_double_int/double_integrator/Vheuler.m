function [t, x, y, w, v, obsx, obsy, obsw, obsv, acc, accy, Ht] = Vheuler(h, t_final, x, y, w, v, obsx, obsy, obsw, obsv, Atol, Rtol, model)
global   vstar n obs_n p  mu2 vmax mu1 A phi a c L lambda   q  rho epsilon visc wmax DYNAMIC_STEP


% t_vec=0:h:t_final;
% F = zeros(n,1);
% u = zeros(n,1);
% K = zeros(n,1);
t = 0;
x_Heun = x; y_Heun = y; w_Heun = w; v_Heun = v;
obs_x = obsx; obs_y = obsy; obs_w = obsw; obs_v = obsv;
% j = 1;

dis = zeros(n);
Vdot = zeros(n); 
t(1) = 0;
z = 1;
vk1(n) = 0; vk2(n) = 0;
wk1(n) = 0; wk2(n) = 0;
xk1(n) = 0; xk2(n) = 0;
yk1(n) = 0; yk2(n) = 0;

% add for obstacles
obs_vk1(obs_n) = 0; obs_vk2(obs_n) = 0;
obs_wk1(obs_n) = 0; obs_wk2(obs_n) = 0;
obs_xk1(obs_n) = 0; obs_xk2(obs_n) = 0;
obs_yk1(obs_n) = 0; obs_yk2(obs_n) = 0;

H(1)=h; 
it = 0;
while (t(z) <= t_final && it < 10000)
    err = 0; error = 0; 
    vk1 = 0; vk2 = 0; wk1 = 0; wk2 = 0; xk1 = 0; xk2 = 0; yk1 = 0; yk2 = 0;
    obs_vk1 = 0; obs_vk2 = 0; obs_wk1 = 0; obs_wk2 = 0; obs_xk1 = 0; obs_xk2 = 0; obs_yk1 = 0; obs_yk2 = 0;
    for i = 1 : n
        V1(i) = 0; V2(i) = 0; viscx(i) = 0; viscy(i) = 0;
        for j = 1 : obs_n
            % if i ~= j
            % dis(i,j) = dist(x_Heun(z,i), x_Heun(z,j), y_Heun(z,i), y_Heun(z,j), p);
            dis(i,j) = dist(x_Heun(z,i), obs_x(z,j), y_Heun(z,i), obs_y(z,j), p);
            Vdot(i,j) = Vd(dis(i,j),  L,lambda, q);
            V1(i) = V1(i) + Vdot(i,j)*(x_Heun(z,i) - obs_x(z,j))/dis(i,j);
            V2(i) = V2(i) + Vdot(i,j)*(y_Heun(z,i) - obs_y(z,j))/dis(i,j);
            if visc ~= 0
                viscx(i) = viscx(i) + kappa(dis(i,j),L,lambda,visc) * (g1(obs_v(z,j)) - g1(v_Heun(z,i)));
                viscy(i) = viscy(i) + kappa(dis(i,j),L,lambda,visc) * (g2(obs_w(z,j)) - g2(w_Heun(z,i)));
            end
            % end
        end
        %%
        if model == 1
            Lambdai = -V1(i) + viscx(i);
            Gi =  Ud(y_Heun(z,i), a, c) + p*V2(i) - viscy(i);
            Li = mu2 + 1/(mu2*wmax)*Gi^2;
            k = mu1 + (Lambdai /vstar)+(vmax / (vstar * (vmax - vstar))) * f(-Lambdai, epsilon);
            F = - k * (v_Heun(z,i) - vstar) + Lambdai;
            u = - Li*w_Heun(z,i) - Gi;
        else
            qC = (vmax*v_Heun(z,i)*cos(w_Heun(z,i))+vstar*vmax-2*vstar*v_Heun(z,i))/(2*(vmax-v_Heun(z,i))^2*v_Heun(z,i)^2);
            betaC = A/((cos(w_Heun(z,i))-cos(phi))^2)+((rho-1)*v_Heun(z,i)*cos(w_Heun(z,i))+vstar)/(vmax-v_Heun(z,i));
            aC = rho*vmax*sin(w_Heun(z,i))/(2*(vmax-v_Heun(z,i))^2*v_Heun(z,i));
            F  = 1/qC*(viscx(i)-mu2*(v_Heun(z,i)*cos(w_Heun(z,i))-vstar)-V1(i));
            u  =  v_Heun(z,i)/betaC*(viscy(i)-mu1*v_Heun(z,i)*sin(w_Heun(z,i))-Ud(y_Heun(z,i),a,c)- aC*F-p*V2(i));
        end
        %%
                        
        acc(z,i) = F;
        accy(z,i) = u;
        
        vk1(i) = h*(F);
        wk1(i) = h*(u);
        yk1(i) = h*(w_Heun(z,i));
        xk1(i) = h*(v_Heun(z,i));
        
        v_temp(i) = v_Heun(z,i) + vk1(i);
        w_temp(i) = w_Heun(z,i) + wk1(i);
        y_temp(i) = y_Heun(z,i) + yk1(i);
        x_temp(i) = x_Heun(z,i) + xk1(i);
        
        % compute obstacle temp trajectories
        for j = 1 : obs_n
            F_obs = 0.0; u_obs = 0.0;
            obs_vk1(j) = h*(F_obs);
            obs_wk1(j) = h*(u_obs);
            obs_yk1(j) = h*(obs_w(z,j));
            obs_xk1(j) = h*(obs_v(z,j));
        
            obs_v_temp(j) = obs_v(z,j) + obs_vk1(j);
            obs_w_temp(j) = obs_w(z,j) + obs_wk1(j);
            obs_y_temp(j) = obs_y(z,j) + obs_yk1(j);
            obs_x_temp(j) = obs_x(z,j) + obs_xk1(j);
        end
    end
    
    for i = 1 : n
        
        V1(i)=0; V2(i)=0;
        viscx(i)=0; viscy(i)=0;
        for j = 1 : obs_n
            dis(i,j) = dist(x_temp(i) , obs_x(z,j) + obs_xk1(j), y_temp(i), obs_y(z,j) + obs_yk1(j), p);
            Vdot(i,j) = Vd(dis(i,j), L, lambda,  q);
            V1(i) = V1(i) + Vdot(i,j)*(x_temp(i) - obs_x(z,j) - obs_xk1(j))/dis(i,j);
            V2(i) = V2(i) + Vdot(i,j)*(y_temp(i) - obs_y(z,j) - obs_yk1(j))/dis(i,j);  
            if visc ~=0
                viscx(i)=viscx(i)+kappa(dis(i,j),L,lambda,visc)*(g1(obs_v(z,j)+vk1(j))-g1(v_Heun(z,i)+vk1(i)));
                viscy(i)=viscy(i)+kappa(dis(i,j),L,lambda,visc)*(g2(obs_w(z,j)+wk1(j))-g2(w_Heun(z,i)+wk1(i)));
            end
        end
        %%
        if model == 1
            Lambdai = -V1(i) + viscx(i);
            Gi =  Ud(y_Heun(z,i) + yk1(i), a, c) + p*V2(i) - viscy(i);
            Li = mu2 + 1 / (mu2*wmax) * Gi^2;
            k = mu1 + (Lambdai /vstar) + (vmax / (vstar * (vmax - vstar))) * f(-Lambdai, epsilon);
            F = - k * (v_Heun(z,i) + vk1(i) - vstar) + Lambdai;
            u = - Li * (w_Heun(z,i) + wk1(i)) - Gi;
        else
            qC = (vmax*(v_Heun(z,i)+vk1(i))*cos(w_Heun(z,i)+wk1(i))+vstar*vmax-2*vstar*(v_Heun(z,i)+vk1(i)))/(2*(vmax-(v_Heun(z,i)+vk1(i)))^2*(v_Heun(z,i)+vk1(i))^2);
            betaC = A/((cos(w_Heun(z,i)+wk1(i))-cos(phi))^2)+((rho-1)*(v_Heun(z,i)+vk1(i))*cos(w_Heun(z,i)+wk1(i))+vstar)/(vmax-(v_Heun(z,i)+vk1(i)));
            aC = rho*vmax*sin(w_Heun(z,i)+wk1(i))/(2*(vmax-(v_Heun(z,i)+vk1(i)))^2*(v_Heun(z,i)+vk1(i)));
            F  = 1/qC*(viscx(i)-mu2*((v_Heun(z,i)+vk1(i))*cos(w_Heun(z,i)+wk1(i))-vstar)-V1(i));
            u  =  (v_Heun(z,i)+vk1(i))/betaC*(viscy(i)-mu1*(v_Heun(z,i)+vk1(i))*sin(w_Heun(z,i)+wk1(i))-Ud(y_Heun(z,i)+yk1(i),a,c)- aC*F-p*V2(i));
        end
        %%
        
        vk2(i) = h*(F);
        wk2(i) = h*(u);
        yk2(i) = h*(w_temp(i));
        xk2(i) = h*(v_temp(i));
        

        %%
        v_Heun(z+1,i) = v_Heun(z,i) + (0.5)*(vk1(i) + vk2(i));
        w_Heun(z+1,i) = w_Heun(z,i) + (0.5)*(wk1(i) + wk2(i));
        y_Heun(z+1,i) = y_Heun(z,i) + (0.5)*(yk1(i) + yk2(i));
        x_Heun(z+1,i) = x_Heun(z,i) + (0.5)*(xk1(i) + xk2(i));

        % 
        for j = 1 : obs_n
            F_obs = 0.0; u_obs = 0.0;
            obs_vk2(j) = h*(F_obs);
            obs_wk2(j) = h*(u_obs);
            obs_yk2(j) = h*(obs_w_temp(j));
            obs_xk2(j) = h*(obs_v_temp(j));

            obs_v(z+1,j) = obs_v(z,j) + (0.5)*(obs_vk1(j) + obs_vk2(j));
            obs_w(z+1,j) = obs_w(z,j) + (0.5)*(obs_wk1(j) + obs_wk2(j));
            obs_y(z+1,j) = obs_y(z,j) + (0.5)*(obs_yk1(j) + obs_yk2(j));
            obs_x(z+1,j) = obs_x(z,j) + (0.5)*(obs_xk1(j) + obs_xk2(j));
        end
         
        %%
        sc_x(i) = Atol + Rtol * max(abs(x_Heun(z+1,i)), abs(x_Heun(z,i)));
        sc_y(i) = Atol + Rtol * max(abs(y_Heun(z+1,i)), abs(y_Heun(z,i)));
        sc_w(i) = Atol + Rtol * max(abs(w_Heun(z+1,i)), abs(w_Heun(z,i)));
        sc_v(i) = Atol + Rtol * max(abs(v_Heun(z+1,i)), abs(v_Heun(z,i)));
        
        err(i) = sqrt(0.25*((x_temp(i) - x_Heun(z+1,i))/sc_x(i))^2 + 0.25*((y_temp(i) - y_Heun(z+1,i))/sc_y(i))^2 + 0.25*((v_temp(i) - v_Heun(z+1,i))/sc_v(i))^2 + 0.25*((w_temp(i) - w_Heun(z+1,i))/sc_w(i))^2);
        %-----------------------------%
        if (i == n && DYNAMIC_STEP == 1)
            error = max(err());
            if error <= 1
                z = z + 1;
                h = h*min(2, 0.8*sqrt(1/error));
                t(z) = t(z-1) + h;
            elseif error > 1
                h = h*min(2, 0.8*sqrt(1/error));
            end
        else
            z = z + 1;
            t(z) = t(z-1) + h;
        end
    end
  
    Ht(z) = h;
    
    it = it + 1;
end
x = transpose(x_Heun); w = transpose(w_Heun); y = transpose(y_Heun); 
v = transpose(v_Heun); acc = transpose(acc); accy = transpose(accy);

obsx = transpose(obs_x); obsw = transpose(obs_w); obsy = transpose(obs_y); obsv = transpose(obs_v);
end

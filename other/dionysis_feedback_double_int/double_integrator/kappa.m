function val = kappa(d,L,lambda,visc)
 
%  
if (L <= d) && (d <= lambda)
    val = visc*(lambda-d)^2;
else
 val = 0;


% if (L <= d) && (d <= lambda)
%     val = (5.225*1e-04)*(lambda-d)^2/(d-L);
% else
%  val = 0;
end
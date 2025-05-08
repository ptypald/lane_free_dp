function val =  Ud(y,a,c)

% val =1*(2*y)/(a^2 - y^2)^2;

% if (y>=-c && y<=c)  WRONG
%     val =0;
% else
%     val=log((a^2-(y-sign(y)*c)^2)/a^2)^2;
% end
% 
if (y>=-a*sqrt(c-1)/sqrt(c) && y<=a*sqrt(c-1)/sqrt(c))
    val = 0;
else
    val = 1*4*(2*y)/(a^2-y^2)^2*(1/(a^2-y^2)-c/a^2)^3;
end  

end
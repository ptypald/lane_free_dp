function val = U(y,a,c)

 
% 
% if (y>=-c && y<=c) WRONG
%     val =0;
% elseif (-a<=y && y<=c) 
%     val=log((a^2-(y -c)^2)/a^2)^2;
% else
%      val=log((a^2-(y +c)^2)/a^2)^2;
% end


 %val =(1)/(a^2 - y^2)^2-1/a^2;
if (y >= -a * sqrt(c-1) / sqrt(c) && y <= a*sqrt(c-1)/sqrt(c))
    val = 0;
else
    val =  (1/(a^2-y^2)-c/a^2)^4;
end

end


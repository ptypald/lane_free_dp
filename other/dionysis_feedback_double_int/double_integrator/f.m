function val =  f(x,epsilon)

 
%  val = epsilon/2+1/(2*epsilon)*x^2;

if  x <= -epsilon
    val=0;
elseif (-epsilon < x && x <= 0)
    val= 1/(2*epsilon)*(x+epsilon)^2;
else
    val = 1/(2*epsilon)*(epsilon^2+2*epsilon*x);
end
%     
    
    
end
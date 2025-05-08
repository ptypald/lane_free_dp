function val = Vd(d,L,lambda,q)
 

if (L <= d) && (d <= lambda)
    val = -q*(d-lambda)^2*(2*d-3*L+lambda)/(d-L)^2;
else
 val = 0;
 
end 
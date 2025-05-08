function val = V(d,L,lambda,q)
 

if (L <= d) && (d <= lambda)
    val = (q)*(lambda-d)^3/(d-L);
else
 val = 0;
end
function val = g1(s)
 
global vmax vstar
 
val = 1/2*((vstar-s)/(s*(s-vmax))+((log((s*(vstar-vmax))/(vstar*(s-vmax))))/(vmax)));
end
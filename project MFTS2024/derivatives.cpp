djxdvx = 0.5*penalty*(1. - ux)*(-LONK3 - (LONK3*(-UXMIN - LONK3*vx))/Sqrt(epsilon + Power(-UXMIN - LONK3*vx,2)))*(ux*UXMAX - uxPrev + 0.5*(1. - ux)*(UXMIN - LONK3*vx + Sqrt(epsilon + Power(-UXMIN - LONK3*vx,2))));

djxdux = 1.*penalty*(UXMAX - 0.5*(UXMIN - LONK3*vx + Sqrt(epsilon + Power(-UXMIN - LONK3*vx,2))))*(ux*UXMAX - uxPrev + 0.5*(1. - ux)*(UXMIN - LONK3*vx + Sqrt(epsilon + Power(-UXMIN - LONK3*vx,2))));



djydy = 1.*penalty*(-(LATK1*(1. - uy)) - LATK1*uy)*(-uyPrev + (1. - uy)*(-(LATK2*(vy - vyLB)) - LATK1*(y - yLB)) + uy*(-(LATK2*(vy - vyUB)) - LATK1*(y - yUB)));

djydvy = 1.*penalty*(-(LATK2*(1. - uy)) - LATK2*uy)*(-uyPrev + (1. - uy)*(-(LATK2*(vy - vyLB)) - LATK1*(y - yLB)) + uy*(-(LATK2*(vy - vyUB)) - LATK1*(y - yUB)));

djyduy = 1.*penalty*(-uyPrev + (1. - uy)*(-(LATK2*(vy - vyLB)) - LATK1*(y - yLB)) + uy*(-(LATK2*(vy - vyUB)) - LATK1*(y - yUB)))*(LATK2*(vy - vyLB) - LATK2*(vy - vyUB) + LATK1*(y - yLB) - LATK1*(y - yUB));



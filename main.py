from hairyBH import black_hole as bh
from hairyBH import particula_time_like as ptl

x=bh(1,12.52,1.4)
t=ptl(1,12.52,1.2,0,2.6e-6)
print(x.masa(), t.U_potencial(0.5))
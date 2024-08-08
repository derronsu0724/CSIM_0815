* Sample HSPICE Netlist
.SUBCKT main n1 n2 n3 n4
C2 n1 0 1
R2 n2 0 2
L2 n3 0 33
L1 n4 0 4
.ENDS main  
V1 1 0 DC 5V
C2 2 1 10n
R1 1 2 1k
X1 1 2 13 4 main
C1 2 0 10n
L2 2 0 11
.tran 0.1ms 10ms
.end

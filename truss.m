% Nodal coordinates
nA = [0; 0];
nB = [5; 4];
nC = [10; 4];
nD = [15; 4];
nE = [20; 0];
nF = [15; 0];
nG = [10; 0];
nH = [5; 0];

% Element Length
L1 = sqrt((nB(1)-nA(1))^2+(nB(2)-nA(2))^2); 
L2 = sqrt((nC(1)-nB(1))^2+(nB(2)-nC(2))^2); 
L3 = sqrt((nD(1)-nC(1))^2+(nD(2)-nC(2))^2); 
L4 = sqrt((nE(1)-nD(1))^2+(nE(2)-nD(2))^2); 
L5 = sqrt((nF(1)-nE(1))^2+(nF(2)-nE(2))^2); 
L6 = sqrt((nG(1)-nF(1))^2+(nG(2)-nF(2))^2); 
L7 = sqrt((nH(1)-nG(1))^2+(nH(2)-nG(2))^2); 
L8 = sqrt((nA(1)-nH(1))^2+(nA(2)-nH(2))^2); 
L9 = sqrt((nB(1)-nH(1))^2+(nB(2)-nH(2))^2); 
L10 = sqrt((nH(1)-nC(1))^2+(nH(2)-nC(2))^2);
L11 = sqrt((nF(1)-nD(1))^2+(nF(2)-nD(2))^2);
L12 = sqrt((nF(1)-nC(1))^2+(nF(2)-nC(2))^2);
L13 = sqrt((nG(1)-nC(1))^2+(nG(2)-nC(2))^2);

% Angle of the elements
th1 = atan((nB(2)-nA(2))/(nB(1)-nA(1)));
th2 = atan((nC(2)-nB(2))/(nC(1)-nB(1)));
th3 = atan((nD(2)-nC(2))/(nD(1)-nC(1)));
th4 = atan((nE(2)-nD(2))/(nE(1)-nD(1)));
th5 = atan((nF(2)-nE(2))/(nF(1)-nE(1)));
th6 = atan((nG(2)-nF(2))/(nG(1)-nF(1)));
th7 = atan((nH(2)-nG(2))/(nH(1)-nG(1)));
th8 = atan((nA(2)-nH(2))/(nA(1)-nH(1)));
th9 = atan((nB(2)-nH(2))/(nB(1)-nH(1)));
th10 = atan((nH(2)-nC(2))/(nH(1)-nC(1)));
th11 = atan((nF(2)-nD(2))/(nF(1)-nD(1)));
th12 = atan((nF(2)-nC(2))/(nF(1)-nC(1)));
th13 = atan((nG(2)-nC(2))/(nG(1)-nC(1)));
% Sines and cosines
C1 = cos(th1);
C2 = cos(th2);
C3 = cos(th3);
C4 = cos(th4);
C5 = cos(th5);
C6 = cos(th6);
C7 = cos(th7);
C8 = cos(th8);
C9 = cos(th9);
C10 = cos(th10);
C11 = cos(th11);
C12 = cos(th12);
C13 = cos(th13);

S1 = sin(th1);
S2 = sin(th2);
S3 = sin(th3);
S4 = sin(th4);
S5 = sin(th5);
S6 = sin(th6);
S7 = sin(th7);
S8 = sin(th8);
S9 = sin(th9);
S10 = sin(th10);
S11 = sin(th11);
S12 = sin(th12);
S13 = sin(th13);
% Transformation matrix
T1 = [C1 S1 0 0; 0 0 C1 S1];
T2 = [C2 S2 0 0; 0 0 C2 S2];
T3 = [C3 S3 0 0; 0 0 C3 S3];
T4 = [C4 S4 0 0; 0 0 C4 S4];
T5 = [C5 S5 0 0; 0 0 C5 S5];
T6 = [C6 S6 0 0; 0 0 C6 S6];
T7 = [C7 S7 0 0; 0 0 C7 S7];
T8 = [C8 S8 0 0; 0 0 C8 S8];
T9 = [C9 S9 0 0; 0 0 C9 S9];
T10 = [C10 S10 0 0; 0 0 C10 S10];
T11 = [C11 S11 0 0; 0 0 C11 S11];
T12 = [C12 S12 0 0; 0 0 C12 S12];
T13 = [C13 S13 0 0; 0 0 C13 S13];

% Cross section of the elements (m^2)
A = 1e-3;
% Material
E = 200e9;
% Boundary conditions
V1 = -30000;
V2 = -60000;
V3 = -30000;

n = 16; % degrees of freedom

% Elemental stiffness matrix in the local coordinate system
% Elelent stiffess
k1 = E*A/L1;
k2 = E*A/L2;
k3 = E*A/L3;
k4 = E*A/L4;
k5 = E*A/L5;
k6 = E*A/L6;
k7 = E*A/L7;
k8 = E*A/L8;
k9 = E*A/L9;
k10 = E*A/L10;
k11 = E*A/L11;
k12 = E*A/L12;
k13 = E*A/L13;
% Stiffness matrix
KL1 = [k1 -k1; -k1 k1];
KL2 = [k2 -k2; -k2 k2];
KL3 = [k3 -k3; -k3 k3];
KL4 = [k4 -k4; -k4 k4];
KL5 = [k5 -k5; -k5 k5];
KL6 = [k6 -k6; -k6 k6];
KL7 = [k7 -k7; -k7 k7];
KL8 = [k8 -k8; -k8 k8];
KL9 = [k9 -k9; -k9 k9];
KL10 = [k10 -k10; -k10 k10];
KL11 = [k11 -k11; -k11 k11];
KL12 = [k12 -k12; -k12 k12];
KL13 = [k13 -k13; -k13 k13];

% Elemental stiffness matrix in the global coordinate system
KG1 = T1'*KL1*T1;
KG2 = T2'*KL2*T2;
KG3 = T3'*KL3*T3;
KG4 = T4'*KL4*T4;
KG5 = T5'*KL5*T5;
KG6 = T6'*KL6*T6;
KG7 = T7'*KL7*T7;
KG8 = T8'*KL8*T8;
KG9 = T9'*KL9*T9;
KG10 = T10'*KL10*T10;
KG11 = T11'*KL11*T11;
KG12 = T12'*KL12*T12;
KG13 = T13'*KL13*T13;

% Elemental stiffness matrix in the global coordinate system (degrees of freedom of the structure)
K1 = zeros(n,n); %AB
K1(1:2,1:2) = KG1(1:2,1:2); 
K1(1:2,3:4) = KG1(1:2,3:4);
K1(3:4,1:2) = KG1(3:4,1:2);
K1(3:4,3:4) = KG1(3:4,3:4);
K2 = zeros(n,n); %BC
K2(3:4,3:4) = KG2(1:2,1:2);
K2(3:4,5:6) = KG2(1:2,3:4);
K2(5:6,3:4) = KG2(3:4,1:2);
K2(5:6,5:6) = KG2(3:4,3:4);
K3 = zeros(n,n); %CD
K3(5:6,5:6) = KG3(1:2,1:2);
K3(5:6,7:8) = KG3(1:2,3:4);
K3(7:8,5:6) = KG3(3:4,1:2);
K3(7:8,7:8) = KG3(3:4,3:4);
K4 = zeros(n,n); %DE
K4(7:8,7:8) = KG4(1:2,1:2);
K4(7:8,9:10) = KG4(1:2,3:4);
K4(9:10,7:8) = KG4(3:4,1:2);
K4(9:10,9:10) = KG4(3:4,3:4);
K5 = zeros(n,n); %EF
K5(9:10,9:10) = KG5(1:2,1:2);
K5(9:10,11:12) = KG5(1:2,3:4);
K5(11:12,9:10) = KG5(3:4,1:2);
K5(11:12,11:12) = KG5(3:4,3:4);
K6 = zeros(n,n); %FG
K6(11:12,11:12) = KG6(1:2,1:2);
K6(11:12,13:14) = KG6(1:2,3:4);
K6(13:14,11:12) = KG6(3:4,1:2);
K6(13:14,13:14) = KG6(3:4,3:4);
K7 = zeros(n,n); %GH
K7(13:14,13:14) = KG7(1:2,1:2);
K7(13:14,15:16) = KG7(1:2,3:4);
K7(15:16,13:14) = KG7(3:4,1:2);
K7(15:16,15:16) = KG7(3:4,3:4);
K8 = zeros(n,n); %HA
K8(15:16,15:16) = KG8(1:2,1:2);
K8(15:16,1:2) = KG8(1:2,3:4);
K8(1:2,15:16) = KG8(3:4,1:2);
K8(1:2,1:2) = KG8(3:4,3:4);
K9 = zeros(n,n); %HB
K9(15:16,15:16) = KG9(1:2,1:2);
K9(15:16,3:4) = KG9(1:2,3:4);
K9(3:4,15:16) = KG9(3:4,1:2);
K9(3:4,3:4) = KG9(3:4,3:4);
K10 = zeros(n,n); %HC
K10(15:16,15:16) = KG10(1:2,1:2);
K10(15:16,5:6) = KG10(1:2,3:4);
K10(5:6,15:16) = KG10(3:4,1:2);
K10(5:6,5:6) = KG10(3:4,3:4);
K11 = zeros(n,n); %FD
K11(11:12,11:12) = KG11(1:2,1:2);
K11(11:12,7:8) = KG11(1:2,3:4);
K11(7:8,11:12) = KG11(3:4,1:2);
K11(7:8,7:8) = KG11(3:4,3:4);
K12 = zeros(n,n); %FC
K12(11:12,11:12) = KG12(1:2,1:2);
K12(11:12,5:6) = KG12(1:2,3:4);
K12(5:6,11:12) = KG12(3:4,1:2);
K12(5:6,5:6) = KG12(3:4,3:4);
K13 = zeros(n,n); %GC
K13(13:14,13:14) = KG13(1:2,1:2);
K13(13:14,5:6) = KG13(1:2,3:4);
K13(5:6,13:14) = KG13(3:4,1:2);
K13(5:6,5:6) = KG13(3:4,3:4);

% Assembly of the stiffness matrix of the structure
K = K1+K2+K3+K4+K5+K6+K7+K8+K9+K10+K11+K12+K13;

% Boundary conditions (displacements)
BCD = zeros(n,n);
BCD(1,1) = 1;
BCD(2,2) = 1;
BCD(10,10) = 1;
D0 = zeros(n,1);

% Boundary conditions (forces)
BCR = zeros(n,n);
BCR(3,3)= 1;
BCR(4,4)= 1;
BCR(5,5)= 1;
BCR(6,6)= 1;
BCR(7,7)= 1;
BCR(8,8)= 1;
BCR(9,9)= 1;
%BCR(10,9)= 1;
BCR(11,11)= 1;
BCR(12,12)= 1;
BCR(13,13)= 1;
BCR(14,14)= 1;
BCR(15,15)= 1;
BCR(16,16)= 1;
R0 = [0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; V3; 0; V2; 0; V1];

% Identity and zero matrix
I = eye(n);
O = zeros(n,1);

% Solution of the equation of the structure
x = [K -I; BCD BCR]\[O ; D0+R0];

% Internal forces (N)
R1 = KL1*T1*[x(1); x(2); x(3); x(4)]
R2 = KL2*T2*[x(3); x(4); x(5); x(6)]
R3 = KL3*T3*[x(5); x(6); x(7); x(8)]
R4 = KL4*T4*[x(7); x(8); x(9); x(10)]
R5 = KL5*T5*[x(9); x(10); x(11); x(12)]
R6 = KL6*T6*[x(11); x(12); x(13); x(14)]
R7 = KL7*T7*[x(13); x(14); x(15); x(16)]
R8 = KL8*T8*[x(15); x(16); x(17); x(18)]
R9 = KL9*T9*[x(15); x(16); x(3); x(4)]
R10 = KL10*T10*[x(15); x(16); x(5); x(6)]
R11 = KL11*T11*[x(11); x(12); x(7); x(8)]
R12 = KL12*T12*[x(11); x(12); x(5); x(6)]
R13 = KL13*T13*[x(13); x(14); x(5); x(6)]
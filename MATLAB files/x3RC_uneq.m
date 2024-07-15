%Finding the objective function
%When we have a 2 sector lattice would unequal capacitances help

n = 1;      %No. of stages
R = 1;    %Total resistance
C = 1;  %Each capacitor
T_d = 2*R/2*C;

R6 = @(r1,r2,r3,r4,r5) 1 - r1 - r2 - r3 - r4 - r5;

f_min = 0.0001; 
f_max = 1.1/8;

num_points_per_decade = 250;
num_decades = log10(f_max) - log10(f_min);
num_points = num_points_per_decade * num_decades;
frequencies = linspace(f_min, f_max, 1000);

w = 2*pi*frequencies;
s = 1i*w;

r_mid = R/4;
func2 = @(x) sinc(x).*exp(-1j*2*pi.*x);


%i6/i9
tfi6 = @(r8, r9, c3, x) (4*r9 + r8 + (1./(1j*2*pi*c3.*x)))./((1./(1j*2*pi*c3.*x)) - r8);
%V3/i9
tfv3 = @(r7, r8, r9, c3, x) 2*r9 + r8 + (2*r7 + r8)*tfi6(r8, r9, c3, x);
%i3/i9
tfi3 = @(r5, r6, r7, r8, r9, c2, c3 , x) (2.*tfv3(r7, r8, r9, c3, x) + tfi6(r8, r9, c3, x).*( 4*r6 + r5 + (1./(1j*2*pi*c2.*x)) ))./((1./(1j*2*pi*c2.*x)) - r5);
%V2/i9
tfv2 = @(r4, r5, r6, r7, r8, r9, c2, c3 , x) tfv3(r7, r8, r9, c3, x) + (2*r6 + r5).*tfi6(r8, r9, c3, x) + (2*r4 + r5).*tfi3(r5, r6, r7, r8, r9, c2, c3 , x);
%i1/i9
tfi1 = @(r2, r3, r4, r5, r6, r7, r8, r9, c, c2, c3 , x) (2.*tfv2(r4, r5, r6, r7, r8, r9, c2, c3 , x) + tfi3(r5, r6, r7, r8, r9, c2, c3 , x).*( 4*r3 + r2 + (1./(1j*2*pi*c.*x)) ))./((1./(1j*2*pi*c.*x)) - r2);

%Final transfer function
func1 = @(r1, r2, r3, r4, r5, r6, r7, r8, r9, c, c2, c3, x) 2./( tfv2(r4, r5, r6, r7, r8, r9, c2, c3 , x) + (2*r3 + r2)*tfi3(r5, r6, r7, r8, r9, c2, c3 , x) + (2*r1 + r2)*tfi1(r2, r3, r4, r5, r6, r7, r8, r9, c, c2, c3 , x) );



objective = @(u) integral(@(x) (abs(func1(u(1), u(2), u(3), u(4), u(5), u(6), u(7), u(8), u(9), u(10), u(11), u(10), x) - func2(x)).^2),0,1.1/8);
objective2 = @(u) max(abs(func1(u(1), u(2), u(3), u(4), u(5), u(6), u(7), u(8), u(9), u(10), u(11), u(12), frequencies) - func2(frequencies)));
objective3 = @(u) max(abs(func1(u(1), u(2), u(1), u(1), u(2), u(1), u(1), u(2), u(1), u(3), u(3), u(3), frequencies) - func2(frequencies)));


options = optimoptions('fminunc', 'Display', 'iter', 'Algorithm', 'quasi-newton');
r1_0 = 1/12;
r2_0 = 1/6;
r3_0 = 1/12;
r4_0 = 1/12;
r5_0 = 1/6;
r6_0 = 1/12;
r7_0 = 1/12;
r8_0 = 1/6;
r9_0 = 1/12;
c1_0 = 0.85;
c2_0 = 0.05;
c3_0 = 0.85;

x1 = [0.0243    0.1797    0.1394    0.1394    0.1745    0.1234    0.1234    0.0869    0.0090    0.9792  0.1755    2.2250];
%4.4298e-07
x2 = [0.0239    0.1799    0.1394    0.1394    0.1747    0.1233    0.1233    0.0872    0.0088    0.9786  0.1750    2.2233];
%4.3855e-07
x3 = [0.2 0.3 0.7 1 0.4 1 0.2 0.04 0.8 1.5 0.5 1];
x0 = [r1_0 r2_0 r3_0 r4_0 r5_0 r6_0 r7_0 r8_0 r9_0 c1_0 c2_0 c3_0];
A = [1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0];
b = 1;
lb = [0.0001 0.0001 0.0001 0.0001 0.0001 0.0001 0.0001 0.0001 0.0001 0.0001 0.0001 0.0001];
ub = [1 1 1 1 1 1 1 1 1];

x01 = [r1_0 r2_0 r3_0 r4_0 r5_0 r6_0 r7_0 r8_0 r9_0 c1_0 c2_0 c3_0];
A1 = [1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0];
b1 = 1;
lb1 = [0.0001 0.0001 0.0001 0.0001 0.0001 0.0001 0.0001 0.0001 0.0001 0.0001 0.0001];
x02 = [0.0156    0.1336    0.1321    0.1321    0.1733    0.1322    0.1322    0.1334    0.0155    1.3888 0.1692 1.3888];
ub1 = [1 1 1 1 1 1 1 1 1 2];


opti2 = fmincon(objective2, x02, [],[],A1,b1,lb1,[]);
opti1 = fmincon(objective2, x02, [],[],A1,b1,lb1,[],[],options);
opti3 = fmincon(objective3, [0 0 1], [],[],[6 3 0],1,[0.0001 0.0001 0.0001],[1 1 2]);
disp(opti1)
disp(opti2)
disp(opti3)

disp(objective2(opti1))
disp(objective2(opti2))
disp(objective3(opti3))
%disp(objective([1.1 1.1 1.1]))


for i = 1:1000
    %arr41(i) = angle(func1(opti1(1),opti1(2),opti1(3),opti1(4),opti1(5),opti1(6),opti1(7),opti1(8),opti1(9),opti1(10),opti1(11),opti1(12),frequencies(i)))*180/pi;
    arrx(i) = angle(func1(x1(1),x1(2),x1(3),x1(4),x1(5),x1(6),x1(7),x1(8),x1(9),x1(10),x1(11),x1(12),frequencies(i)))*180/pi;
    %arr430(i) = abs(func1(opti1(1),opti1(2),opti1(3),opti1(4),opti1(5),opti1(6),opti1(7),opti1(8),opti1(9),opti1(10),opti1(11),opti1(12),frequencies(i)) - func2(frequencies(i)));
    arr43(i) = abs(func1(x1(1),x1(2),x1(3),x1(4),x1(5),x1(6),x1(7),x1(8),x1(9),x1(10),x1(11),x1(12),frequencies(i)) - func2(frequencies(i)));
end

hold on
%plot(frequencies, arr43, LineWidth=3)
hold off
legend("Optimized",fontsize=20)
xlabel('Frequency (in Hz)',FontSize=20);
ylabel('Magnitude of Difference',FontSize=20);
%ylabel('Phase (in degrees)',FontSize=20);
%title('Difference in magnitude for optimum',FontSize=20)
%title('Phase Response',FontSize=20)
title('Difference in magnitude for optimum',FontSize=20)
ax = gca;
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
grid on
%writematrix(arr41, 'uneq31.dat')
disp(objective(opti1))

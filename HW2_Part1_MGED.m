clear
clc
close all

%% Input data
load=readtable('load.csv');
loads=load.NYDemand_MWh_(1:168);
scaled_load=(loads-mean(loads))/std(loads)+5;
weather=readtable("weather.csv");
wind_speed=weather.windspeed(1:168);
temperature=weather.temp(1:168);

if weather.solarradiation(1:168) == 0
   solar=weather.solarradiation(1:168);
else
   solar=weather.solarradiation(1:168)+70;
end
%% Wind 1.5 MW
[row,column]=size(loads);
v_ms = wind_speed; %wind speed initialize
%fills v_ms with values from that table 
for i = 1:row%row to 100
    v_ms(1,i) = wind_speed(i);
end
v_ms=v_ms+3;  %scale the wind speed
v_i = 5;%cut in
v_0 = 45;%cut out
v_nom = 15;%nominal wind speed
p_nom1 = 1.5;%MW
p_nom2 = 3;%MW

Pw1 = zeros(1, row);
Pw2 = zeros(1, row);

for i = 1:row%
    if(v_ms(i) < v_i) || (v_ms(1,i) > v_0)
        Pw1(i) = 0;
        Pw2(i) = 0;
    elseif (v_i<= v_ms(i)) && (v_ms(1, i) <= v_nom)
        Pw1(1,i) = p_nom1*((v_ms(1,i)-v_i)/(v_nom - v_i));
        Pw2(1,i) = p_nom2*((v_ms(1,i)-v_i)/(v_nom - v_i));
    elseif (v_nom <= v_ms(1,i)) && (v_ms(i) < v_0)
        Pw1(1,i) = p_nom1;
        Pw2(1,i) = p_nom2;
    end
end

%% Solar 1 & 2 MW
for i=1:168
Pvnom1 = 1000; %1 MW
Pvnom2 = 2000; %1 MW
Pv1(i)=(Pvnom1*solar(i)*(1-0.0047*(temperature(i)-25))/1000)/1000;
Pv2(i)=(Pvnom2*solar(i)*(1-0.0047*(temperature(i)-25))/1000)/1000;
end

%% Cost Coefficients 
% combined heat and power (CHP)
%parameter of a, b, c  
a_chp = 1530;
b_chp = 0.01;
c_chp = 0.000233;
maxChp = 6;
RChp = 0.6;

% RChp = 1;
% diesel generation
a_ds = 1300;
b_ds = 0.013;
c_ds = 0.000235;
maxDs = 4;
RDs = 0.5;

% natural gas
a_ng = 992;
b_ng = 0.016;
c_ng = 0.000241;
maxNg = 10;
RNg = 1;

%cost parameter of Pv, Wind, Batt
k_pv   = 0.0018;
k_wind = 0.0025;
k_bat=0.025;

% cost of 
c_up=0.001;
c_down=0.001;
%% Batteries 
Pb1max = 4; %MW
Pb2max = 6; %MWW
Capb1 = 8; %MWH
Capb2 = 12; %MWH
effb = 0.95;
E0b1 = Capb1;
E0b2 = Capb2;
Pb10 = 0;
Pb20 = 0;
%% Pgrid
Pgbound = 1000; %MW
Cbuy = [20; 20; 25; 27.5; 29; 36; 36.5; 37; 40; 42; 45; 46; 42; 40; 37.5; 36.5; 36; 32; 27.5; 27; 23; 20.5; 20; 20]; %hourly data c/kWh for 1 day
priceBase = 0.005;
Cbuyd = priceBase*Cbuy*1000; %$/Mwh 
daily_Csell = [30.00 29.75 27.00 27.50 30.00 31.00 36.00]; %$/Mwh
combdata=[Cbuyd' daily_Csell];
Cbuyd_m=mean(Cbuyd);
mul=Cbuyd/Cbuyd_m;
Csel=zeros(24,7);
for i =1:7
 Csel(:,i)=daily_Csell(i).*mul;
end
% Final 
Cselw=reshape(Csel,1,168)/1000;
Cbuyw= Cselw/0.8;
%% Initialization
Tfinal = 168;
ts = 1;
nvar = 21;
x0 = [zeros(18,1); E0b1; E0b2; 0];
Pchp = x0(1);
Uchp = x0(2);
SdnChp = x0(3);
SupChp = x0(4);
Pds = x0(5);
Uds = x0(6);
SdnDs = x0(7);
SupDs = x0(8);
Png = x0(9);
Ung = x0(10); 
SdnNg = x0(11);
SupNg = x0(12);
Pw10 = x0(13);
Pw20 = x0(14);
Pv10 = x0(15);
Pv20 = x0(16);
Pb1 = x0(17);
Pb2 = x0(18);
E0b1 = x0(19);
E0b2 = x0(20);
Pgrid = x0(21);
results = zeros(168,22); %1st col Pload, 2nd col Pg1 and so on until Pwind
results(:,1) = scaled_load;
%% Optimization
for i = 1:Tfinal

% x1    x2   x3    x4      x5 x6    x7    x8  x9  x10   x11   x12  x13 x14 x15x16 x17 x18 x19 x20 x21
% Pchp Uchp SdnChp SupChp Pds Uds SdnDs SupDs Png Ung SdnNg SupNg Pw1 Pw2 Pv1 Pv2 Pb1 Pb2 Eb1 Eb2 Pgrid 
        % x1    x2   x3    x4      x5 x6    x7    x8  x9  x10   x11   x12  x13 x14 x15x16   x17   x18       x19 x20 x21
        % Pchp Uchp SdnChp SupChp Pds Uds SdnDs SupDs Png Ung SdnNg SupNg Pw1 Pw2 Pv1 Pv2   Pb1   Pb2       Eb1 Eb2 Pgrid 
Aeq =      [1    0     0      0    1   0    0      0   1    0   0       0  1    1   1   1    1      1       0   0   1; % balance
            0    0     0      0    0   0    0      0   0    0   0       0  0    0   0   0   effb    0       1   0   0;%Energy bat 1
            0    0     0      0    0   0    0      0   0    0   0       0  0    0   0   0     0     effb    0   1   0];%energy bat 2
beq = [scaled_load(i); E0b1; E0b2];


    % x1    x2      x3    x4      x5    x6     x7    x8    x9    x10   x11   x12  x13 x14 x15x16   x17       x18       x19 x20 x21
    % Pchp Uchp  SdnChp SupChp Pds   Uds     SdnDs SupDs Png    Ung    SdnNg SupNg Pw1 Pw2 Pv1 Pv2   Pb1    Pb2       Eb1 Eb2 Pgrid 
A = [-1   -RChp     0      0      0    0       0      0    0      0     0       0  0    0   0   0     0       0       0     0    0; %chp rampdown
     1      0       0      0      0    0       0      0    0      0     0       0  0    0   0   0     0       0       0     0    0;%chp rampup
     0      0       0      0      -1  -RDs     0      0    0      0     0       0  0    0   0   0     0       0       0     0    0; %ds rampdown
     0      0       0      0      1    0       0      0    0      0     0       0  0    0   0   0     0       0       0     0    0;%ds rampup
     0      0       0      0      0    0       0      0    -1   -RNg    0       0  0    0   0   0     0       0       0     0    0;%ng rampdown
     0      0       0      0      0    0       0      0     1      0    0       0  0    0   0   0     0       0       0     0    0;%ng rampup
     -1     1*0     0      0      0    0       0      0    0      0     0       0  0    0   0   0     0       0       0     0    0;%chp lower bound
     1      -maxChp 0      0      0    0       0      0    0      0     0       0  0    0   0   0     0       0       0     0    0;%chp upper bound
     0      0       0      0      -1   1*0     0      0    0      0     0       0  0    0   0   0     0       0       0     0    0;%ds lower bound
     0      0       0      0      1   -maxDs   0      0    0      0     0       0  0    0   0   0     0       0       0     0    0;%ds upper bound
     0      0       0      0      0    0       0      0    -1   1*0     0       0  0    0   0   0     0       0       0     0    0;%ng lower bound
     0      0       0      0      0    0       0      0    1    -maxNg  0       0  0    0   0   0     0       0       0     0    0];%ng upper bound
 b = [-Pchp;Pchp+Uchp*RChp;-Pds;Pds+Uds*RDs;-Png;Png+Ung*RNg;0;0;0;0;0;0];

 %nonlinear constraint
 nlcon = @(x) [x(4) - max(0,x(2)-Uchp);
               x(3) - max(0,-x(2)+Uchp);
               x(8) - max(0,x(6)-Uds);
               x(7) - max(0,-x(6)+Uds);
               x(12) - max(0,x(10)-Ung);
               x(11) - max(0,-x(10)+Ung);];
nlrhs = [0;0;0;0;0;0];
nle = [0;0;0;0;0;0];
%Boundries
% x1    x2   x3    x4      x5 x6    x7    x8  x9  x10   x11   x12  x13 x14 x15`x16 x17 x18  x19 x20 x21
% Pchp Uchp SdnChp SupChp Pds Uds SdnDs SupDs Png Ung SdnNg SupNg Pw1 Pw2 Pv1 Pv2 Pb1 Pb2   Eb1 Eb2 Pgrid 
lb = [zeros(16,1); -Pb1max; -Pb2max; 0.15*Capb1; 0.15*Capb2; 0];
ub = [maxChp; 1; 1; 1; maxDs; 1; 1; 1; maxNg; 1; 1; 1; Pw1(i); Pw2(i); Pv1(i); Pv2(i); Pb1max; Pb2max; Capb1; Capb2; 1000];

%integer Constraints
xtype = 'CBBBCBBBCBBBCCCCCCCCC';

%objfun
fun = @(x) c_chp*x(1)^2 + b_chp*x(1) + a_chp ...
    + c_ds*x(5)^2  + b_ds*x(5)  + a_ds...
    + c_ng*x(9)^2  + b_ng*x(9)  + a_ng ...
    + k_pv*x(15) +k_pv*x(16) + k_wind*x(13)+k_wind*x(14) + k_bat*x(17)+k_bat*x(18)...
    +c_up*x(4)+c_down*x(3)+c_up*x(8)+c_down*x(7)+c_up*x(12)+c_down*x(11)...
    +(Cbuyw(i)+Cselw(i))*x(21)/2+(Cbuyw(i)-Cselw(i))*abs(x(21))/2;

% Create OPTI Object
Opt = opti('fun',fun,'nlmix',nlcon,nlrhs,nle,'ineq',A,b,'eq',Aeq,beq,'bounds',lb,ub,...
    'xtype',xtype);

% Solve the MINLP problem
[x,fval,exitflag,info] = solve(Opt,x0);

results(i,2:22) = x';

%check power balance
pbal(i,1) = beq(1,1)-x(1)-x(5)-x(9)-x(13)-x(14)-x(15)-x(16)-x(17)-x(18)-x(21);
SOCb1(i) = ((E0b1 - effb*x(17))/Capb1)*100;
SOCb2(i) = ((E0b2 - effb*x(18))/Capb2)*100;

%update values in constraints 
Pchp = x(1);
Uchp = x(2);
SdnChp = x(3);
SupChp = x(4);
Pds = x(5);
Uds = x(6);
SdnDs = x(7);
SupDs = x(8);
Png = x(9);
Ung = x(10); 
SdnNg = x(11);
SupNg = x(12);
Pw10 = x(13);
Pw20 = x(14);
Pv10 = x(15);
Pv20 = x(16);
Pb1 = x(17);
Pb2 = x(18);
E0b1 = x(19);
E0b2 = x(20);
Pgrid = x(21);
end

%% Input Data
figure, 
subplot(2,2,1)
plot(Pv1,'y-*','LineWidth',1,'MarkerSize',4)
hold on
plot(Pv2,'-s','LineWidth',1,'Color','#ffa500','MarkerSize',4)
% title('PV Power','fontsize',12,'interpreter','latex')
legend('$P_{PV1}$','$P_{PV2}$','Interpreter','latex','FontSize',12, 'Location','best','Orientation','horizontal')
xlim([1 168])
ylim([-0.5 2.5])
xlabel('Time (hour)','fontsize',12,'interpreter','latex')
ylabel('$P_{PVi}$ (MW)','fontsize',12,'interpreter','latex')
grid on
ax = gca;
ax.GridAlpha = 0.3;
ax.GridLineStyle = '--';

subplot(2,2,2)
stairs(Pw1,'g-^','LineWidth',1,'MarkerSize',4)
hold on
stairs(Pw2,'-o','LineWidth',1,'Color','#279461','MarkerSize',4)
xlim([1 168])
ylim([-0.5 3.5])
% title('Wind Power','fontsize',12,'interpreter','latex')
legend('$P_{W1}$','$P_{W2}$','Interpreter','latex','FontSize',12, 'Location','best','Orientation','horizontal')
xlabel('Time (hour)','fontsize',12,'interpreter','latex')
ylabel('$P_{Wind} (MW)$','fontsize',12,'interpreter','latex')
grid on
ax = gca;
ax.GridAlpha = 0.3;
ax.GridLineStyle = '--';

subplot(2,2,3)
plot(Cbuyw,'m-^','LineWidth',1,'MarkerSize',4)
hold on
plot(Cselw,'c-o','LineWidth',1,'MarkerSize',4)
xlim([1 168])
% title('Dynamic price','fontsize',12,'interpreter','latex')
legend('$C_{buy}$','$C_{sell}$','Interpreter','latex','FontSize',12, 'Location','best','Orientation','horizontal')
xlabel('Time (hour)','fontsize',12,'interpreter','latex')
ylabel('$Cost (\$/kWh)$','fontsize',12,'interpreter','latex')
grid on
ax = gca;
ax.GridAlpha = 0.3;
ax.GridLineStyle = '--';

subplot(2,2,4)
plot(scaled_load,'r-.o','LineWidth',1,'MarkerSize',4)
% hold on
% bar(noise(1,1:length(pload)),'LineStyle','--','EdgeColor','cyan','FaceColor','blue')
title('Load Profile','fontsize',12,'interpreter','latex')
xlim([1 168])
% legend('$P_{load}$','noise','Interpreter','latex','FontSize',12)
xlabel('Time (hour)','fontsize',12,'interpreter','latex')
ylabel('$P_{load}$ (MW)','fontsize',12,'interpreter','latex')
grid on
ax = gca;
ax.GridAlpha = 0.3;
ax.GridLineStyle = '--';
%%

figure,
subplot(3,2,1)
allpd = [results(:,2) results(:,6) results(:,10) results(:,14) results(:,15) results(:,16) results(:,17) results(:,18) results(:,19) results(:,22)];
% bar(allpd,'stacked')
area(allpd)
hold on
plot(results(:,1),'r-d','LineWidth',2,'MarkerSize',4)
title('Power generation','fontsize',12,'interpreter','latex','FontSize',12)
xlabel('Time (hour)','fontsize',12,'interpreter','latex','FontSize',12)
xlim([1 168])
ylabel('$P_i$ (MW)','fontsize',12,'interpreter','latex','FontSize',12)
legend('$P_{CHP}$','$P_{DS}$','$P_{NG}$','$P_{W1}$','$P_{W2}$','$P_{PV1}$','$P_{PV2}$','$P_{B1}$','$P_{B2}$','$P_{Grid}$','$P_{Load}$'...
    ,'Interpreter','latex','FontSize',12, 'Location','best','Orientation','vertical')
grid on
ax = gca;
ax.GridAlpha = 1;
ax.GridLineStyle = ':';

subplot(3,2,5)
stairs(SOCb1,'-o','LineWidth',1,'Color','#00BFBF','MarkerSize',4)
hold on
stairs(SOCb2,'-s','LineWidth',1,'Color','#BF00BF','MarkerSize',4)
title('State of charge of batteries','fontsize',12,'interpreter','latex')
xlabel('Time (hour)','fontsize',12,'interpreter','latex','FontSize',12)
xlim([1 168])
ylabel('$SOC$ ($\%$)','fontsize',12,'interpreter','latex','FontSize',12)
legend('$SOC_{B1}$','$SOC_{B2}$','Interpreter','latex','FontSize',12)
grid on
ax = gca;
ax.GridAlpha = 1;
ax.GridLineStyle = ':';

subplot(3,2,2)
plot(results(:,2),'--','LineWidth',1,'Color','#ff7f01','Marker','v','MarkerSize',4)
hold on
plot(results(:,6),'--s','LineWidth',1,'Color','#142B8C','Marker','hexagram','MarkerSize',4)
hold on
plot(results(:,10),'-.o','LineWidth',1,'Color','#00BFBF','MarkerSize',4)
hold on
plot(results(:,22),'r-.','LineWidth',1,'MarkerSize',4)
title('Output power of generators','fontsize',12,'interpreter','latex')
xlabel('Time (hour)','fontsize',12,'interpreter','latex','FontSize',12)
xlim([1 168])
ylabel('$P_{Gi}$ (MW)','fontsize',12,'interpreter','latex','FontSize',12)
legend('$P_{CHP}$','$P_{DS}$','$P_{NG}$','$P_{Grid}$','Interpreter','latex','FontSize',12, 'Location','best','Orientation','horizontal','FontSize',12)
grid on
ax = gca;
ax.GridAlpha = 1;
ax.GridLineStyle = ':';

subplot(3,2,3)
plot(results(:,14),'g-d','LineWidth',1)
hold on
plot(results(:,15),'-s','LineWidth',1,'Color','#279461')
title('Wind turbines output power','fontsize',12,'interpreter','latex')
xlabel('Time (hour)','fontsize',12,'interpreter','latex')
ylabel('$P_{Wi}$ (MW)','fontsize',12,'interpreter','latex')
xlim([1 168])
legend('$P_{W1}$','$P_{W2}$','Interpreter','latex','FontSize',12, 'Location','best','Orientation','horizontal')
grid on
ax = gca;
ax.GridAlpha = 1;
ax.GridLineStyle = ':';

subplot(3,2,4)
plot(results(:,16),'y-*','LineWidth',1)
hold on
plot(results(:,17),'-s','LineWidth',1,'Color','#ffa500')
title('Solar PV output power','fontsize',12,'interpreter','latex','FontSize',12)
xlabel('Time (hour)','fontsize',12,'interpreter','latex','FontSize',12)
xlim([1 168])
ylabel('$P_{PVi}$ (MW)','fontsize',12,'interpreter','latex','FontSize',12)
legend('$P_{PV1}$','$P_{PV2}$','Interpreter','latex','FontSize',12, 'Location','best','Orientation','horizontal','FontSize',12)
grid on
ax = gca;
ax.GridAlpha = 1;
ax.GridLineStyle = ':';

subplot(3,2,6)
plot(results(:,18),'-o','LineWidth',1,'Color','#00BFBF','MarkerSize',4)
hold on
plot(results(:,19),'-s','LineWidth',1,'Color','#BF00BF','MarkerSize',4)
title('BESS output power','fontsize',12,'interpreter','latex','FontSize',12)
xlabel('Time (hour)','fontsize',12,'interpreter','latex','FontSize',12)
xlim([1 168])
ylabel('$P_{Bi}$ (MW)','fontsize',12,'interpreter','latex','FontSize',12)
legend('$P_{B1}$','$P_{B2}$','Interpreter','latex','FontSize',12, 'Location','best','Orientation','horizontal','FontSize',12)
grid on
ax = gca;
ax.GridAlpha = 1;
ax.GridLineStyle = ':';

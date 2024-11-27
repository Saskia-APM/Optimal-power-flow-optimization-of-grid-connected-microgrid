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

Pw1 = zeros(1, row); %holds values for new wind over month row to 100
Pw2 = zeros(1, row);

for i = 1:row %goes through the list that is the size of v_ms row to 100
%     wind_power(1,i) = 4;
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
Pv1(1,80:92) = Pv1(1,152:164);
Pv2(1,80:92) = Pv2(1,152:164);

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
daily_Csell = [77.91 66.04 64.56 68.00 63.00 74.10 94.65]; %$/Mwh
combdata=[Cbuyd' daily_Csell];
Cbuyd_m=mean(Cbuyd);
mul=Cbuyd/Cbuyd_m;
Csel=zeros(24,7);
for i =1:7
 Csel(:,i)=daily_Csell(i).*mul;
end
Cselw=reshape(Csel,1,168);

%% 
Y =  [(1/1i)+(1/1i) -(1/1i)     -(1/1i);
      -(1/1i)     (1/1i)+(1/1i) -(1/1i);
      -(1/1i)      -(1/1i)    (1/1i)+(1/1i)];
theta = angle(Y);
deg = theta*57.296;

%% Initialization
Tfinal = 168;
ts = 1;
nvar = 35;
x0 = [zeros(18,1); E0b1; E0b2;zeros(15,1)];
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
del1 = x0(21);
del2 = x0(22);
del3 = x0(23);
V1 = x0(24);
V2 = x0(25);
V3 = x0(26);
Qchp=x0(27);
Qds=x0(28);
Qng=x0(29);
Qw1=x0(30);
Qw2=x0(31);
Qpv1=x0(32);
Qpv2=x0(33);
Qb1=x0(34);
Qb2=x0(35);

results = zeros(168,37); %1st col Pload, 2nd col Pg1 and so on until Pwind
results(:,1) = scaled_load;
qload = scaled_load*0.2;
results(:,37) = qload;
%% Optimization
for i = 1:Tfinal

            % x1    x2   x3    x4      x5 x6    x7    x8  x9  x10   x11   x12  x13 x14 x15x16 x17 x18 x19 x20 x21
            % Pchp Uchp SdnChp SupChp Pds Uds SdnDs SupDs Png Ung SdnNg SupNg Pw1 Pw2 Pv1 Pv2 Pb1 Pb2 Eb1 Eb2 Pgrid 
Aeq =      [0    0     0      0    0   0    0      0   0    0   0       0  0    0   0   0   effb    0       1   0       0   0       0   0 0  0  0   0  0 0  0  0  0  0  0;%Energy bat 1
            0    0     0      0    0   0    0      0   0    0   0       0  0    0   0   0     0     effb    0   1       0   0       0   0 0  0  0   0  0 0  0  0  0  0  0];%energy bat 2
beq = [E0b1; E0b2];


%ramp up rampdown constraints 
    % x1    x2      x3    x4      x5    x6     x7    x8    x9    x10   x11 x12  x13 x14 x15x16   x17        x18      x19 x20    x21 x22 x23 x24 x25 x26 x27 x28 x29 x30 x31 x32
    % Pchp Uchp  SdnChp SupChp Pds   Uds     SdnDs SupDs Png    Ung SdnNg SupNg Pw1 Pw2 Pv1 Pv2   Pb1         Pb2    Eb1 Eb2    PN1 PN2 PN3 QN1 QN2 QN3 del1 del2 del3 V1 V2 V3 
A = [-1   -RChp     0      0      0    0       0      0    0      0     0       0  0    0   0   0     0       0       0     0    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; %chp rampdown
     1    -RChp     0      0      0    0       0      0    0      0     0       0  0    0   0   0     0       0       0     0    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;%chp rampup
     0      0       0      0      -1  -RDs     0      0    0      0     0       0  0    0   0   0     0       0       0     0    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; %ds rampdown
     0      0       0      0       1  -RDs     0      0    0      0     0       0  0    0   0   0     0       0       0     0    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;%ds rampup
     0      0       0      0      0    0       0      0    -1   -RNg    0       0  0    0   0   0     0       0       0     0    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;%ng rampdown
     0      0       0      0      0    0       0      0     1   -RNg    0       0  0    0   0   0     0       0       0     0    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;%ng rampup
     -1     1*0     0      0      0    0       0      0    0      0     0       0  0    0   0   0     0       0       0     0    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;%chp lower bound
     1      -maxChp 0      0      0    0       0      0    0      0     0       0  0    0   0   0     0       0       0     0    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;%chp upper bound
     0      0       0      0      -1   1*0     0      0    0      0     0       0  0    0   0   0     0       0       0     0    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;%ds lower bound
     0      0       0      0      1   -maxDs   0      0    0      0     0       0  0    0   0   0     0       0       0     0    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;%ds upper bound
     0      0       0      0      0    0       0      0    -1   1*0     0       0  0    0   0   0     0       0       0     0    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;%ng lower bound
     0      0       0      0      0    0       0      0    1    -maxNg  0       0  0    0   0   0     0       0       0     0    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%ng upper bound
 b = [-Pchp;Pchp;-Pds;Pds;-Png;Png;0;0;0;0;0;0];
%           x1   x2   x3    x4     x5 x6    x7    x8  x9  x10   x11  x12  x13 x14 x15x16   x17   x18       x19 x20     x21
        % Pchp Uchp SdnChp SupChp Pds Uds SdnDs SupDs Png Ung SdnNg SupNg Pw1 Pw2 Pv1 Pv2   Pb1   Pb2       Eb1 Eb2     Pg1 Pg2 Pg3
 %nonlinear constraint
 nlcon = @(x) [x(4) - max(0,x(2)-Uchp);
               x(3) - max(0,-x(2)+Uchp);
               x(8) - max(0,x(6)-Uds);
               x(7) - max(0,-x(6)+Uds);
               x(12) - max(0,x(10)-Ung);
               x(11) - max(0,-x(10)+Ung);
               x(1)+x(15)+x(13)-abs(x(24))*abs(x(24))*abs(Y(1,1))*cos(theta(1,1))-abs(x(24))*abs(x(25))*abs(Y(1,2))*cos(theta(1,2)-x(21)+x(22))-abs(x(24))*abs(x(26))*abs(Y(1,3))*cos(theta(1,3)-x(21)+x(23));
               x(9)+x(16)+x(17)-scaled_load(i)/2-abs(x(25))*abs(x(24))*abs(Y(2,1))*cos(theta(2,1)-x(22)+x(21))-abs(x(25))*abs(x(25))*abs(Y(2,2))*cos(theta(2,2))-abs(x(25))*abs(x(26))*abs(Y(2,3))*cos(theta(2,3)-x(22)+x(23));
               x(5)+x(14)+x(18)-scaled_load(i)/2-abs(x(26))*abs(x(24))*abs(Y(3,1))*cos(theta(3,1)-x(23)+x(21))-abs(x(26))*abs(x(25))*abs(Y(3,2))*cos(theta(3,2)-x(23)+x(22))-abs(x(26))*abs(x(26))*abs(Y(3,3))*cos(theta(3,3));
               x(27)+x(30)+x(32)-abs(x(24))*abs(x(27))*abs(Y(1,1))*sin(-theta(1,1))-abs(x(24))*abs(x(25))*abs(Y(1,2))*sin(-theta(1,2)+x(21)-x(22))-abs(x(24))*abs(x(26))*abs(Y(1,3))*sin(-theta(1,3)+x(21)-x(23));
               x(29)+x(33)+x(24)-(scaled_load(i)*0.2/2)-abs(x(25))*abs(x(24))*abs(Y(2,1))*sin(-theta(2,1)+x(22)-x(21))-abs(x(25))*abs(x(25))*abs(Y(2,2))*sin(-theta(2,2))-abs(x(24))*abs(x(26))*abs(Y(2,3))*sin(-theta(2,3)+x(22)-x(23));
               x(28)+x(35)+x(31)-(scaled_load(i)*0.2/2)-abs(x(26))*abs(x(24))*abs(Y(3,1))*sin(-theta(3,1)+x(23)-x(21))-abs(x(26))*abs(x(25))*abs(Y(3,2))*sin(-theta(3,2)+x(23)-x(22))-abs(x(26))*abs(x(26))*abs(Y(3,3))*sin(-theta(3,3))];

nlrhs = [0;0;0;0;0;0;0;0;0;0;0;0];
nle = [0;0;0;0;0;0;0;0;0;0;0;0];

%Boundaries
% x1    x2   x3    x4      x5 x6    x7    x8  x9  x10   x11   x12  x13 x14 x15`x16 x17 x18  x19 x20 x21
% Pchp Uchp SdnChp SupChp Pds Uds SdnDs SupDs Png Ung SdnNg SupNg Pw1 Pw2 Pv1 Pv2 Pb1 Pb2   Eb1 Eb2 Pgrid 
lb = [zeros(16,1); -Pb1max; -Pb2max; 0.2*Capb1; 0.2*Capb2;-pi;-pi;-pi;5.4;5.4;5.4;-maxChp;-maxDs;-maxNg;-Pw1(i);-Pw2(i);-Pv1(i);-Pv2(i);-Pb1max;-Pb2max];
ub = [maxChp; 1; 1; 1; maxDs; 1; 1; 1; maxNg; 1; 1; 1; Pw1(i); Pw2(i); Pv1(i); Pv2(i); Pb1max; Pb2max; Capb1; Capb2;...
    pi;pi;pi;6.6;6.6;6.6;maxChp;maxDs;maxNg;Pw1(i);Pw2(i);Pv1(i);Pv2(i);Pb1max;Pb2max];

%type of dec 
xtype = 'CBBBCBBBCBBBCCCCCCCCCCCCCCCCCCCCCCC';

%objfun
fun = @(x) c_chp*x(1)^2 + b_chp*x(1) + a_chp ...
    + c_ds*x(5)^2  + b_ds*x(5)  + a_ds...
    + c_ng*x(9)^2  + b_ng*x(9)  + a_ng ...
    + k_pv*x(15) +k_pv*x(16) + k_wind*x(13)+k_wind*x(14) + k_bat*x(17)+k_bat*x(18)...
    +c_up*x(4)+c_down*x(3)+c_up*x(8)+c_down*x(7)+c_up*x(12)+c_down*x(11);


% Create OPTI Object
Opt = opti('fun',fun,'nlmix',nlcon,nlrhs,nle,'ineq',A,b,'eq',Aeq,beq,'bounds',lb,ub,...
    'xtype',xtype);


% Solve the MINLP problem
[x,fval,exitflag,info] = solve(Opt,x0);


results(i,2:36) = x';

 %check power balance
pbal(i,1) = beq(1,1)-x(1)-x(5)-x(9)-x(13)-x(14)-x(15)-x(16)-x(17)-x(18);
pg1bal(i,1) = x(1) + x(15) + x(13)-scaled_load(i)/3;
pg2bal(i,1) = x(9) + x(16) + x(17)-scaled_load(i)/3;
pg3bal(i,1) = x(5) + x(18) + x(14)-scaled_load(i)/3;
pb11= pg1bal(:,1)+ pg2bal(:,1)+pg3bal(:,1);

EF(i,1) = exitflag;
SOCb1(i) = ((E0b1 - effb*x(17))/Capb1)*100;
SOCb2(i) = ((E0b2 - effb*x(18))/Capb2)*100;

 %update values
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
del1 = x(21);
del2 = x(22);
del3 = x(23);
V1 = x(24);
V2 = x(25);
V3 = x(26);
Qchp=x(27);
Qds=x(28);
Qng=x(29);
Qw1=x(30);
Qw2=x(31);
Qpv1=x(32);
Qpv2=x(33);
Qb1=x(34);
Qb2=x(35);
end
%% data configuration
PG1 = results(:,2) + results(:,16) + results(:,14);
PG2 = results(:,10) + results(:,17) + results(:,18);
PG3 = results(:,6) + results(:,15) + results(:,19);
QG1 = results(:,28) + results(:,31) + results(:,33);
QG2 = results(:,30) + results(:,34) + results(:,25);
QG3 = results(:,29) + results(:,36) + results(:,32);
resultsQ = [results(:,28)  results(:,31)  results(:,33) results(:,30)  results(:,34)  results(:,25) results(:,29)  results(:,36)  results(:,32)];
%% Input data
figure('name','name','Units','inches',...
  'Position',[1 1 5 3],...
  'PaperPositionMode','auto');
plot(Pv1,'y-*','LineWidth',2,'Color','#ffd700','MarkerIndices',1:2:length(Pv1))
hold on
plot(Pv2,'r-x','LineWidth',2,'Color','#ffa500','MarkerIndices',1:2:length(Pv1))
% title('PV Power','fontsize',12,'interpreter','latex')
legend('$P_{PV1}$','$P_{PV2}$','Interpreter','latex','FontSize',12, 'Location','best','Orientation','horizontal')
xlim([1 168])
xlabel('Time (hour)','fontsize',12,'interpreter','latex')
ylabel('$P_{PVi}$ (MW)','fontsize',12,'interpreter','latex')
grid on
ax = gca;
ax.GridAlpha = 0.3;
ax.GridLineStyle = '--';
exportgraphics(gca,'PVpower.png','Resolution',600)
%%
figure('name','name','Units','inches',...
  'Position',[1 1 5 3],...
  'PaperPositionMode','auto');
stairs(Pw1,'g-^','LineWidth',1,'MarkerSize',4)
hold on
stairs(Pw2,'-o','LineWidth',1,'Color','#279461','MarkerSize',4)
xlim([1 168])
% title('Wind Power','fontsize',12,'interpreter','latex')
legend('$P_{W1}$','$P_{W2}$','Interpreter','latex','FontSize',12, 'Location','best','Orientation','horizontal')
xlabel('Time (hour)','fontsize',12,'interpreter','latex')
ylabel('$P_{Wind} (MW)$','fontsize',12,'interpreter','latex')
grid on
ax = gca;
ax.GridAlpha = 1;
ax.GridLineStyle = ':';

%%
figure('name','name','Units','inches',...
  'Position',[1 1 5 3],...
  'PaperPositionMode','auto');
plot(Cbuyw,'m-^','LineWidth',1)
hold on
plot(Cselw,'c-o','LineWidth',1)
xlim([1 168])
% title('Dynamic price','fontsize',12,'interpreter','latex')
legend('$C_{buy}$','$C_{sell}$','Interpreter','latex','FontSize',12, 'Location','best','Orientation','horizontal')
xlabel('Time (hour)','fontsize',12,'interpreter','latex')
ylabel('$Cost (\$/MW)$','fontsize',12,'interpreter','latex')
grid on
ax = gca;
ax.GridAlpha = 1;
ax.GridLineStyle = ':';

subplot(2,2,4)
%%
figure('name','name','Units','inches',...
  'Position',[1 1 5 3],...
  'PaperPositionMode','auto');
plot(scaled_load,'r-o','LineWidth',2,'MarkerIndices',1:2:length(scaled_load))
% hold on
% bar(noise(1,1:length(pload)),'LineStyle','--','EdgeColor','cyan','FaceColor','blue')
% title('Load Profile','fontsize',12,'interpreter','latex')
xlim([1 168])
% legend('$P_{load}$','noise','Interpreter','latex','FontSize',12)
xlabel('Time (hour)','fontsize',14,'interpreter','latex')
ylabel('$P_{load}$ (MW)','fontsize',14,'interpreter','latex')
grid on
ax = gca;
ax.GridAlpha = 1;
ax.GridLineStyle = '--';
% exportgraphics(gca,'loadprof.png','Resolution',600)
%%
figure,
subplot(3,2,1)
allpd = [results(:,2) results(:,6) results(:,10) results(:,14) results(:,15) results(:,16) results(:,17) results(:,18) results(:,19)];
area(allpd)
hold on
plot(results(:,1),'r-d','LineWidth',2)
title('Economic Dispatch','fontsize',12,'interpreter','latex','FontSize',12)
xlabel('Time (hour)','fontsize',12,'interpreter','latex','FontSize',12)
xlim([1 168])
ylabel('$Power$ (MW)','fontsize',12,'interpreter','latex','FontSize',12)
legend('$P_{CHP}$','$P_{DS}$','$P_{NG}$','$P_{W1}$','$P_{W2}$','$P_{PV1}$','$P_{PV2}$','$P_{B1}$','$P_{B2}$','Interpreter','latex','FontSize',12, 'Location','best','Orientation','horizontal')
grid on
ax = gca;
ax.GridAlpha = 1;
ax.GridLineStyle = ':';

subplot(3,2,5)
stairs(SOCb1,'-o','LineWidth',1,'Color','#00BFBF')
hold on
stairs(SOCb2,'-s','LineWidth',1,'Color','#BF00BF')
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
plot(results(:,2),'b-','LineWidth',1,'Marker','v')
hold on
plot(results(:,6),'--s','LineWidth',1,'Color','#D95319','Marker','hexagram')
hold on
plot(results(:,10),'k-.o','LineWidth',1)
title('Output power of generators','fontsize',12,'interpreter','latex')
xlabel('Time (hour)','fontsize',12,'interpreter','latex','FontSize',12)
xlim([1 168])
ylabel('$Power$ (MW)','fontsize',12,'interpreter','latex','FontSize',12)
legend('$P_{CHP}$','$P_{DS}$','$P_{NG}$','Interpreter','latex','FontSize',12, 'Location','best','Orientation','horizontal','FontSize',12)
grid on
ax = gca;
ax.GridAlpha = 1;
ax.GridLineStyle = ':';

subplot(3,2,3)
plot(results(:,14),'g-^','LineWidth',1,'MarkerSize',4)
hold on
plot(results(:,15),'-s','LineWidth',1,'Color','#279461','MarkerSize',4)
title('Wind turbines output power','fontsize',12,'interpreter','latex')
xlabel('Time (hour)','fontsize',12,'interpreter','latex')
ylabel('$Power$ (MW)','fontsize',12,'interpreter','latex')
xlim([1 168])
ylim([-0.5 3.5])
legend('$P_{W1}$','$P_{W2}$','Interpreter','latex','FontSize',12, 'Location','best','Orientation','horizontal')
grid on
ax = gca;
ax.GridAlpha = 1;
ax.GridLineStyle = ':';

subplot(3,2,4)
plot(results(:,16),'y-*','LineWidth',1,'MarkerSize',4)
hold on
plot(results(:,17),'-s','LineWidth',1,'Color','#ffa500','MarkerSize',4)
title('Solar PV output power','fontsize',12,'interpreter','latex','FontSize',12)
xlabel('Time (hour)','fontsize',12,'interpreter','latex','FontSize',12)
xlim([1 168])
ylabel('$Power$ (MW)','fontsize',12,'interpreter','latex','FontSize',12)
legend('$P_{PV1}$','$P_{PV2}$','Interpreter','latex','FontSize',12, 'Location','best','Orientation','horizontal','FontSize',12)
grid on
ax = gca;
ax.GridAlpha = 1;
ax.GridLineStyle = ':';

subplot(3,2,6)
plot(results(:,18),'-o','LineWidth',1,'Color','#00BFBF')
hold on
plot(results(:,19),'-s','LineWidth',1,'Color','#BF00BF')
title('Battery output power','fontsize',12,'interpreter','latex','FontSize',12)
xlabel('Time (hour)','fontsize',12,'interpreter','latex','FontSize',12)
xlim([1 168])
ylabel('$Power$ (MW)','fontsize',12,'interpreter','latex','FontSize',12)
legend('$P_{B1}$','$P_{B2}$','Interpreter','latex','FontSize',12, 'Location','best','Orientation','horizontal','FontSize',12)
grid on
ax = gca;
ax.GridAlpha = 1;
ax.GridLineStyle = ':';

%% daily
figure,
subplot(3,2,1)
allpd = [results(:,2) results(:,6) results(:,10) results(:,14) results(:,15) results(:,16) results(:,17) results(:,18) results(:,19)];

area(allpd)
hold on
plot(results(:,1),'r-d','LineWidth',2,'Marker','hexagram')
title('Economic Dispatch','fontsize',12,'interpreter','latex','FontSize',12)
xlabel('Time (hour)','fontsize',12,'interpreter','latex','FontSize',12)
xlim([145 168])
% xlim([1 24])
ylabel('$Power$ (MW)','fontsize',12,'interpreter','latex','FontSize',12)
legend('$P_{CHP}$','$P_{DS}$','$P_{NG}$','$P_{W1}$','$P_{W2}$','$P_{PV1}$','$P_{PV2}$','$P_{B1}$','$P_{B2}$','Interpreter','latex','FontSize',12, 'Location','best','Orientation','horizontal')
grid on
ax = gca;
ax.GridAlpha = 1;
ax.GridLineStyle = ':';

% Set new x-axis ticks and labels for the subplot
new_xticks = linspace(145, 168, 24);
new_xticklabels = arrayfun(@num2str, 1:24, 'UniformOutput', false);
ax.XTick = new_xticks;
ax.XTickLabel = new_xticklabels;

subplot(3,2,5)
stairs(SOCb1(1,:),'-o','LineWidth',2,'Color','#00BFBF')
hold on
stairs(SOCb2(1,:),'-s','LineWidth',2,'Color','#BF00BF')
title('State of charge of batteries','fontsize',12,'interpreter','latex')
xlabel('Time (hour)','fontsize',12,'interpreter','latex','FontSize',12)
xlim([145 168])
% xlim([1 24])
ylabel('$SOC$ ($\%$)','fontsize',12,'interpreter','latex','FontSize',12)
legend('$SOC_{B1}$','$SOC_{B2}$','Interpreter','latex','FontSize',12)
grid on
ax = gca;
ax.GridAlpha = 1;
ax.GridLineStyle = ':';

% Set new x-axis ticks and labels for the subplot
new_xticks = linspace(120, 144, 24);
new_xticklabels = arrayfun(@num2str, 1:24, 'UniformOutput', false);
ax.XTick = new_xticks;
ax.XTickLabel = new_xticklabels;
% ,'Color','#00BFBF'


subplot(3,2,2)
plot(results(:,2),'b-','LineWidth',2,'Marker','v')
hold on
plot(results(:,6),'--s','LineWidth',2,'Color','#D95319','Marker','hexagram')
hold on
plot(results(:,10),'k-.o','LineWidth',2)
title('Output power of generators','fontsize',12,'interpreter','latex')
xlabel('Time (hour)','fontsize',12,'interpreter','latex','FontSize',12)
xlim([145 168])
% xlim([1 24])
ylabel('$Power$ (MW)','fontsize',12,'interpreter','latex','FontSize',12)
legend('$P_{CHP}$','$P_{DS}$','$P_{NG}$','Interpreter','latex','FontSize',12, 'Location','best','Orientation','horizontal','FontSize',12)
grid on
ax = gca;
ax.GridAlpha = 1;
ax.GridLineStyle = ':';

% Set new x-axis ticks and labels for the subplot
new_xticks = linspace(145, 168, 24);
new_xticklabels = arrayfun(@num2str, 1:24, 'UniformOutput', false);
ax.XTick = new_xticks;
ax.XTickLabel = new_xticklabels;

subplot(3,2,3)
plot(results(:,14),'g-^','LineWidth',2,'MarkerSize',4)
hold on
plot(results(:,15),'-o','LineWidth',2,'Color','#279461','MarkerSize',4)
title('Wind turbines output power','fontsize',12,'interpreter','latex')
xlabel('Time (hour)','fontsize',12,'interpreter','latex')
ylabel('$$P_{Wi}$$ (MW)','fontsize',12,'interpreter','latex')
xlim([145 168])
% xlim([1 24])
legend('$P_{W1}$','$P_{W2}$','Interpreter','latex','FontSize',12, 'Location','best','Orientation','horizontal')
grid on
ax = gca;
ax.GridAlpha = 1;
ax.GridLineStyle = ':';
% Set new x-axis ticks and labels for the subplot
new_xticks = linspace(145, 168, 24);
new_xticklabels = arrayfun(@num2str, 1:24, 'UniformOutput', false);
ax.XTick = new_xticks;
ax.XTickLabel = new_xticklabels;

subplot(3,2,4)
plot(results(:,16),'y-*','LineWidth',2,'MarkerSize',4)
hold on
plot(results(:,17),'-s','LineWidth',2,'Color','#ffa500','MarkerSize',4)
title('Solar PV output power','fontsize',12,'interpreter','latex','FontSize',12)
% xlabel('Time (hour)','fontsize',12,'interpreter','latex','FontSize',12)
xlim([145 168])
% xlim([1 24])
ylabel('$P_{PVi}$ (MW)','fontsize',12,'interpreter','latex','FontSize',12)
legend('$P_{PV1}$','$P_{PV2}$','Interpreter','latex','FontSize',12, 'Location','best','Orientation','horizontal','FontSize',12)
grid on
ax = gca;
ax.GridAlpha = 1;
ax.GridLineStyle = ':';

% Set new x-axis ticks and labels for the subplot
new_xticks = linspace(145, 168, 24);
new_xticklabels = arrayfun(@num2str, 1:24, 'UniformOutput', false);
ax.XTick = new_xticks;
ax.XTickLabel = new_xticklabels;

subplot(3,2,6)
plot(results(:,18),'-o','LineWidth',2,'Color','#00BFBF')
hold on
plot(results(:,19),'-s','LineWidth',2,'Color','#BF00BF')
title('Battery output power','fontsize',12,'interpreter','latex','FontSize',12)
xlabel('Time (hour)','fontsize',12,'interpreter','latex','FontSize',12)
xlim([145 168])
% xlim([1 24])
ylabel('$Power$ (MW)','fontsize',12,'interpreter','latex','FontSize',12)
legend('$P_{B1}$','$P_{B2}$','Interpreter','latex','FontSize',12, 'Location','best','Orientation','horizontal','FontSize',12)
grid on
ax = gca;
ax.GridAlpha = 1;
ax.GridLineStyle = ':';

% % Set new x-axis ticks and labels for the subplot
new_xticks = linspace(145, 168, 24);
new_xticklabels = arrayfun(@num2str, 1:24, 'UniformOutput', false);
ax.XTick = new_xticks;
ax.XTickLabel = new_xticklabels;
%% daily
figure,
subplot(3,1,1)
plot(PG1,'b-.o','LineWidth',2)
hold on
plot(PG2,'g--d','LineWidth',2)
hold on
plot(PG3,'-s','LineWidth',2,'Color','#BF00BF')
xlim([145 168])
% xlim([1 24])
title('Active power','fontsize',12,'interpreter','latex','FontSize',12)
xlabel('Time (hour)','fontsize',12,'interpreter','latex','FontSize',12)
ylabel('$P_i$ (MW)','fontsize',12,'interpreter','latex','FontSize',12)
legend('$P_{G1}$','$P_{G2}$','$P_{G3}$','Interpreter','latex','FontSize',12, 'Location','best','Orientation','horizontal')
grid on
ax = gca;
ax.GridAlpha = 1;
ax.GridLineStyle = ':';
% % Set new x-axis ticks and labels for the subplot
new_xticks = linspace(145, 168, 24);
new_xticklabels = arrayfun(@num2str, 1:24, 'UniformOutput', false);
ax.XTick = new_xticks;
ax.XTickLabel = new_xticklabels;

subplot(3,1,2)
plot(QG1,'b-.o','LineWidth',2)
hold on
plot(QG2,'g--d','LineWidth',2)
hold on
plot(QG3,'-s','LineWidth',2,'Color','#BF00BF')
xlim([145 168])
% xlim([1 24])
title('Reactive power','fontsize',12,'interpreter','latex')
xlabel('Time (hour)','fontsize',12,'interpreter','latex')
ylabel('$Q_i$ (MVar)','fontsize',12,'interpreter','latex')
legend('$Q_{1}$','$Q_{2}$','$Q_{3}$','Interpreter','latex','FontSize',12, 'Location','best','Orientation','horizontal')
grid on
ax = gca;
ax.GridAlpha = 1;
ax.GridLineStyle = ':';

% % Set new x-axis ticks and labels for the subplot
new_xticks = linspace(145, 168, 24);
new_xticklabels = arrayfun(@num2str, 1:24, 'UniformOutput', false);
ax.XTick = new_xticks;
ax.XTickLabel = new_xticklabels;

subplot(3,1,3)
plot(results(:,25),'b-.o','LineWidth',2)
hold on
plot(results(:,26),'g--d','LineWidth',2)
hold on
plot(results(:,27),'-s','LineWidth',2,'Color','#BF00BF')
hold on
plot(6.6*ones(168,1),'k--','LineWidth',1)
hold on
plot(5.4*ones(168,1),'k--','LineWidth',1)
title('Voltage profile','fontsize',12,'interpreter','latex')
xlabel('Time (hour)','fontsize',12,'interpreter','latex')
ylabel('$V_i$ (kV)','fontsize',12,'interpreter','latex')
xlim([145 168])
% xlim([1 24])
legend('$V_{1}$','$V_{2}$','$V_{3}$','Interpreter','latex','FontSize',12, 'Location','best','Orientation','horizontal')
grid on
ax = gca;
ax.GridAlpha = 1;
ax.GridLineStyle = ':';
% % Set new x-axis ticks and labels for the subplot
new_xticks = linspace(145, 168, 24);
new_xticklabels = arrayfun(@num2str, 1:24, 'UniformOutput', false);
ax.XTick = new_xticks;
ax.XTickLabel = new_xticklabels;

%% weekly
figure,
subplot(3,1,1)
plot(PG1,'b-.d','LineWidth',1,'MarkerSize',4)
hold on
plot(PG2,'g--o','LineWidth',1,'MarkerSize',4)
hold on
plot(PG3,'-s','LineWidth',1,'Color','#BF00BF','MarkerSize',4)
% hold on
% plot(results(1:24,1),'r-d','LineWidth',2)
title('Active power','fontsize',12,'interpreter','latex','FontSize',12)
% xlabel('Time (hour)','fontsize',12,'interpreter','latex','FontSize',12)
ylabel('$P_i$ (MW)','fontsize',12,'interpreter','latex','FontSize',12)
legend('$P_{G1}$','$P_{G2}$','$P_{G3}$','Interpreter','latex','FontSize',12, 'Location','best','Orientation','horizontal')
grid on
xlim([1 168])
ax = gca;
ax.GridAlpha = 1;
ax.GridLineStyle = ':';

subplot(3,1,2)
plot(QG1,'b-.d','LineWidth',1,'MarkerSize',4)
hold on
plot(QG2,'g--d','LineWidth',1,'MarkerSize',4)
hold on
plot(QG3,'-s','LineWidth',1,'Color','#BF00BF','MarkerSize',4)
hold on
plot(results(:,37),'k--s','LineWidth',1)
title('Reactive power','fontsize',12,'interpreter','latex')
% xlabel('Time (hour)','fontsize',12,'interpreter','latex')
ylabel('$Q_i$ (MVar)','fontsize',12,'interpreter','latex')
legend('$Q_{1}$','$Q_{2}$','$Q_{3}$','Interpreter','latex','FontSize',12, 'Location','best','Orientation','horizontal')
grid on
xlim([1 168])
ax = gca;
ax.GridAlpha = 1;
ax.GridLineStyle = ':';

subplot(3,1,3)
plot(results(:,25),'b-.d','LineWidth',1,'MarkerSize',4)
hold on
plot(results(:,26),'g--o','LineWidth',1,'MarkerSize',4)
hold on
plot(results(:,27),'-s','LineWidth',1,'Color','#BF00BF','MarkerSize',4)
hold on
plot(6.6*ones(168,1),'k--','LineWidth',1)
hold on
plot(5.4*ones(168,1),'k--','LineWidth',1)
title('Voltage profile','fontsize',12,'interpreter','latex')
xlabel('Time (hour)','fontsize',12,'interpreter','latex')
ylabel('$V_i$ (kV)','fontsize',12,'interpreter','latex')
xlim([1 168])
ylim([5.3 6.7])
legend('$V_{1}$','$V_{2}$','$V_{3}$','Interpreter','latex','FontSize',12, 'Location','best','Orientation','horizontal')
grid on
ax = gca;
ax.GridAlpha = 1;
ax.GridLineStyle = ':';
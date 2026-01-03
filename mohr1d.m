% Ustav konstruovania a strojarskych technologii­
% Technicka fakulta
% SPU Nitra
% Vytvorene v GNU Octave 7.1.0
% Zadanie c.1 Pruznost a pevnost

clc;
clearvars;

%******************************
%Zadane hodnoty - sem sa zadavaju hodnoty pre vypocet
sx=50e6;
sy=0;
tz=0;
a=19.896;
%*******************************
%Vypocet
s=((sx+sy)/2)+(((sx-sy)/2)*cosd(2*a))+tz*sind(2*a);
t=(((sx-sy)/2)*sind(2*a))-tz*cosd(2*a);
s1=(sx+sy)/2+0.5*(sqrt((sx-sy)^2+4*tz^2));
s2=(sx+sy)/2-0.5*(sqrt((sx-sy)^2+4*tz^2));
t1=0.5*(sqrt((sx-sy)^2+4*tz^2));
t2=-0.5*(sqrt((sx-sy)^2+4*tz^2));
sk=(s1+s2)/2; %stred kruznice
r=sqrt((0.5*(sx-sy))^2+tz^2); %polomer kruznice
%Vypis vypocitanych hodnot


disp("Vstupne napatia");
disp("sx=");
disp(sx);
disp("sy=");
disp(sy);
disp("tz=");
disp(tz);
disp("Vypocitane hodnoty");
disp("sigma, tau");
disp(s);
disp(t);
disp("Hlavne napatia s1, s2");
disp(s1);
disp(s2);
disp("napatia t1,t2");
disp(t1);
disp(t2);
disp("Polomer kruznice")
disp(r);

srx=r;
sry=0;

% vypocet hodnot kruznice
tt=linspace(0,2*pi,100);
cx=r*cos(tt)+sk;
cy=r*sin(tt);
figure('Name','Mohrová kružnica','NumberTitle','off');
set(gcf,'color','w');
set(gca,'FontSize',14);
hold on;
grid on;
%Vykreslenie kruznice+vypisovanie popiskov
plot(cx,cy);
axis equal;

if(sy==0) % Priamkova napatost

title("Mohrová kružnica-priamková napätost");
xlabel('\sigma ,Pa','interpreter', 'tex',"fontsize", 16);
ylabel('\tau ,Pa','interpreter', 'tex',"fontsize", 16);

plot(sk,0,"r*"); % stred kruznice
text(sk,sk/10,"S","fontsize", 16);
hold on;
plot(s,t,"r*");
text(s,t,' \sigma , tau','interpreter', 'tex',"fontsize", 16);
%text(sy,-tz/10,' \sigmay ','interpreter', 'tex',"fontsize", 16);
hold on;
line([sk s], [0 t], "linestyle", "-", "color", "b");
hold on;
line([sk sx-s],[0 t],"linestyle", "-.", "color", "b");
hold on;
plot(s1,0,"r*");%Sigma1
text(s1,0,' \sigmax','interpreter', 'tex',"fontsize", 16);


hold on;

else
title("Mohrová kružnica-rovinná napätost");
xlabel('\sigma ,Pa','interpreter', 'tex',"fontsize", 16);
ylabel('\tau ,Pa','interpreter', 'tex',"fontsize", 16);
hold on;
plot(sk,0,"r*"); % stred kruznice
text(sk,sk/10,"S","fontsize", 16);
hold on;
plot(s2,0,"r*");%Sigma2
text(s2,0,' \sigma2','interpreter', 'tex',"fontsize", 16);
hold on;
plot(s1,0,"r*");%Sigma1
text(s1,0,' \sigma1','interpreter', 'tex',"fontsize", 16);
hold on;
plot(s,t,"r*");
text(s,t,' \sigma , tau','interpreter', 'tex',"fontsize", 16);
text(sy,-tz/10,' \sigmay ','interpreter', 'tex',"fontsize", 16);
text(sx,-tz/10,' \sigmax ','interpreter', 'tex',"fontsize", 16);
hold on;
line([sk s], [0 t], "linestyle", "-", "color", "b");
hold on;
%ciara cez stred rovnobezne s tau
x=[sk sk];
y=[t1 t2];
line(x,y, "linestyle", "--", "color", "g");
hold on
%vynasanie tz
x1=[sx sx];
y1=[0 -tz];
line(x1,y1, "linestyle", "--", "color", "r");
text(sx,-tz,' \tauz','interpreter', 'tex',"fontsize", 16);
hold on;
%vynasanie tz
x2=[sy sy];
y2=[0 tz];
line(x2,y2, "linestyle", "--", "color", "r");
text(sy,tz,' \tauz','interpreter', 'tex',"fontsize", 16);
%diagonala smyku
hold on;
x3=[sy sx];
y3=[ tz -tz];
line(x3,y3, "linestyle", "-.", "color", "r");
hold on;
%Pol kruznice
x4=[sy sy];
y4=[ 0 -tz];
line(x4,y4, "linestyle", ":", "color", "r");
x5=[sy sx];
y5=[ -tz -tz];
line(x5,y5, "linestyle", ":", "color", "r");
hold on;
plot(sy,-tz,"r*");% Pol kruznice
text(sy,-(tz-tz/10),'P',"fontsize", 16,"color", "m");
hold on;
%Smery hlavnych napatia
x6=[sy s1];
y6=[ -tz 0];
line(x6,y6, "linestyle", "-", "color", "m");
text(s1,-tz/10,'smer \sigma1','interpreter', 'tex',"fontsize", 16,"color", "m");
x7=[sy s2];
y7=[ -tz 0];
line(x7,y7, "linestyle", "-", "color", "m");
text(s2,-tz/10,'smer \sigma2','interpreter', 'tex',"fontsize", 16,"color", "m");
%Vykreslenia sigma, tau=f(phi)
phi=zeros(1,360);
sxx=zeros(1,360);
tau=zeros(1,360);
ph=-180;
for i=1:360
  if (i==1)
    phi(i)=ph;
  else
    phi(i)=ph+i;
  endif
   sxx(i)=0.5*(sx+sy)+0.5*(sx-sy)*cosd(2*phi(i))+tz*sind(2*phi(i));
   tau(i)=0.5*(sx-sy)*sind(2*phi(i))-tz*cosd(2*phi(i));
endfor;
figure('Name','Priebeh napätí','NumberTitle','off');
set(gcf,'color','w');
set(gca,'FontSize',14);
plot(phi,sxx,"linewidth", 1);
set(gca, "linewidth", 1, "fontsize", 14)
hold on;
plot(phi,tau,"linewidth", 1);
hold on;
grid on;
title("Napätia-rovinná napätost","fontsize", 16);
xticks([-180 -100  0 100 180]);
xticklabels({'-180','-100','0','100','180'});
xlabel('\phi ,Deg','interpreter', 'tex',"fontsize", 16);
ylabel(' \sigma, \tau, Pa','interpreter', 'tex',"fontsize", 16);
hold on;
legend('{\fontsize{12} Normálove napätie }', '{\fontsize{12} Šmykové napätie }');

endif



disp("-Koniec vypoctu-");
disp("Viac info:https://www.tecgraf.puc-rio.br/etools/mohr/mohreng.html#signconv");
hold off;





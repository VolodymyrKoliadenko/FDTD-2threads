% TM !!!!!!!!!!!!!
clc; clear all;
%clo=clock
dx=0.5*1e-7/1.4;
dy=dx;% НЕ МЕНЯЙ, нужно переписывать код
period=1.2e-6;%2
Q=1.5;
n=1;%ne используется в Эпсилон матрице пока =1
lambd=1064e-9;
prodol=2*n*period^2/lambd/Q;
nx=1700;%600
ny=900;%1700
Imeem__Sloy= ny*dy/prodol*2
Ex=zeros(nx,ny)+1i*1e-100;Ex0=Ex;%с ноликом текущие, без нолика предыдущие. Чтоб запутать:)
Ey=zeros(nx,ny)+1i*1e-100;Ey0=Ey;
Hz=zeros(nx-1,ny-1)+1i*1e-100;Hz0=Hz;
e=ones(nx,ny);%epsilon
i=2;j=2;k=2;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5


omega=2*pi/lambd;%частота HERZ /3e8
 
dt=dx/2;%сразу тау=ст было 0,2
tau=50000*dt;
%dz=1;
Z=376.7303;
s=dt/dx;  k3=(1-s)/(1+s);
w=15*1e-7;% NE shagov
ii=1:ny;
%_____      E P S I L O N
 
m=(0.008*4*2);% Эпсилон в 2раза больше должно быть /АМПЛИТУДА КАК ВСЕГДА
nachalo=3;
% ПЕРЕНОРМИРОВКА МОДУЛЯЦИИ
for i=1:nx
    e(i,1:nachalo-1)=1;
    e(i,nachalo:end)=1+(m/2)*(1+sign(-0.1+cos(2*pi*(i-nx/2+0.5)*dx/period)*sin(2*pi*(ii(nachalo:end)-nachalo)*dy/prodol)));
end
e=(dt*Z/dx)./e;
dH=dt/dx/Z;
%e=(dt*Z/dx)./e;
fi=0*asind(lambd/2/period);%gradus

clo=clock;
%i=2: ЦИКЛ и й к
tmax=ny*2.2;
for t=1:tmax
    tt=min(t*s+10,ny-1);
    for i=2:nx-1
        
            Ex(i,1)=exp(-((i-nx/2+0.5)^2*dx^2)/w^2+1i*omega*tand(fi)*(i-nx/2+0.5)*dx)*exp(1i*(t-1)*dt*omega)*exp(-(t-1)*dt/tau);% затухания
            %Hz(i,1)=(1/Z)*exp(-((i-n/2+0.5)^2*dx^2)/w^2+1i*omega*tand(fi)*(i-n/2+0.5)*dx)*exp(1i*(t-0.5)*dt*omega + 1i*omega*dy/2)*exp(-(t-0.5)*dt/tau);% затухания
        
    end
for i=2:nx-1
    for j=2:tt
        
            Ex0(i,j)=Ex(i,j)+e(i,j)*((Hz(i,j)-Hz(i,j-1)));
            Ey0(i,j)=Ey(i,j)-e(i,j)*(Hz(i,j)-Hz(i-1,j));
           
            
        
    end
end
for i=1:nx
    Ex0(i,end)= Ex(i,end-1)+k3*(Ex(i,end)-Ex0(i,end-1)); %GOOD задний край
end
for i=1:ny
    Ex0(1,i)= Ex(2,i)+k3*(Ex(1,i)-Ex0(2,i));
    Ey0(1,i)= Ey(2,i)+k3*(Ey(1,i)-Ey0(2,i));
    Ey0(end,i)=Ey(end-1,i)+k3*(Ey(end,i)-Ey0(end-1,i));
end        

 Ex=Ex0; Ey=Ey0; 
% %АШ будет на 1/2 новее
% %аш использует уже посчитанное Е
% %ЦИКЛ и й к
    for i=2:nx-1
        
            Ex(i,1)=exp(-((i-nx/2+0.5)^2*dx^2)/w^2+1i*omega*tand(fi)*(i-nx/2+0.5)*dx)*exp(1i*(t)*dt*omega)*exp(-(t)*dt/tau);%без затухания
            %Hz(i,1)=(1/Z)*exp(-((i-n/2+0.5)^2*dx^2)/w^2+1i*omega*tand(fi)*(i-n/2+0.5)*dx)*exp(1i*(t-0.5)*dt*omega + 1i*omega*dy/2)*exp(-(t-0.5)*dt/tau);% затухания
        
    end
for i=1:nx-1
    for j=1:tt
        % константы упрощения
            Hz0(i,j)=Hz(i,j)+dH*(Ey(i,j)-Ey(i+1,j)-Ex(i,j)+Ex(i,j+1));
        
    end
end
Hz=Hz0;

  if (mod(t,10)==0)
      subplot(1,2,1)
    imagesc((abs(Ex)).^2)
    colormap('hot')
    %colorbar
    c=clock-clo;
    sec=c(6)+60*(c(6)<0);
    time=c(3)*86400+c(4)*3600+c(5)*60+c(6);
    speed=t/time*10;
    title([' Step ' int2str(t) '. Time ' int2str(c(5)-(c(6)<0)) ':' int2str(sec) ' Speed ' int2str(speed) '/10']);
    subplot(1,2,2)
    spectrTemp=single(abs(fftshift(fft(fftshift(Ex)))));
    imagesc(spectrTemp(nx/2-50:nx/2+50,:))
    %clear spectrTemp;
    colormap('hot')
    pause(0.003);
end
end




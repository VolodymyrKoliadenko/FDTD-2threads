clc; clear all; %close all;
%clo=clock
dx=1*1e-7/1.2;
dy=dx;% НЕ МЕНЯЙ, нужно переписывать код
period=2e-6;
Q=1.08;
n=1;%ne используется в Эпсилон матрице пока =1
lambd=1064e-9;%1064
prodol=2*n*period^2/lambd/Q;
nx=512;%600
ny=1000;%1700
Imeem__Sloy= ny*dy/prodol*2
Ez=zeros(nx,ny)+1i*1e-100;Ez0=Ez;%с ноликом текущие, без нолика предыдущие. Чтоб запутать:)
Hy=zeros(nx-1,ny-1)+1i*1e-100;Hy0=Hy;
Hx=zeros(nx-1,ny-1)+1i*1e-100;Hx0=Hx;
e=ones(nx,ny);%epsilon
i=2;j=2;k=2;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5


omega=2*pi/lambd;%частота HERZ /3e8
 
dt=dx/2;%сразу тау=ст было 0,2
tau=50000*dt;
%dz=1;
Z=376.7303;
s=dt/dx;  k3=(1-s)/(1+s);
w=30*1e-7;% NE shagov
ii=1:ny;
%_____      E P S I L O N
 
m=(0.003*2*2);% Эпсилон в 2раза больше должно быть /АМПЛИТУДА КАК ВСЕГДА
nachalo=3;
% ПЕРЕНОРМИРОВКА МОДУЛЯЦИИ
for i=1:nx
    e(i,1:nachalo-1)=1;
    e(i,nachalo:end)=1+(m/2)*(1+sign(-0.1+cos(2*pi*(i-nx/2+0.5)*dx/period)*sin(2*pi*(ii(nachalo:end)-nachalo)*dy/prodol)));
end
e=(dt*Z/dx)./e;
dH=dt/dx/Z;
%e=(dt*Z/dx)./e;
fi=0;%asind(lambd/2/period);%gradus

clo=clock;
figure
%i=2: ЦИКЛ и й к
tmax=ny*2.2;
for t=1:tmax
    tt=int32(min(t*s+10,ny-1));
    for i=2:nx-1
        
            Ez(i,1)=exp(-((i-nx/2+0.5)^2*dx^2)/w^2+1i*omega*tand(fi)*(i-nx/2+0.5)*dx)*exp(1i*(t-1)*dt*omega)*exp(-(t-1)*dt/tau);% затухания
            %Hx(i,1)=(1/Z)*exp(-((i-n/2+0.5)^2*dx^2)/w^2+1i*omega*tand(fi)*(i-n/2+0.5)*dx)*exp(1i*(t-0.5)*dt*omega + 1i*omega*dy/2)*exp(-(t-0.5)*dt/tau);% затухания
        
    end

    % MAIN EZ
    Ez0(2:end-1, 2:tt) = Ez(2:end-1, 2:tt)+e(2:end-1, 2:tt).*(Hx(2:end, 1:tt-1)-Hx(2:end, 2:tt)+Hy(2:end, 2:tt)-Hy(1:end-1, 2:tt));

for i=1:nx
    Ez0(i,end)= Ez(i,end-1)+k3*(Ez(i,end)-Ez0(i,end-1)); %GOOD задний край
end
for i=1:ny
    Ez0(1,i)= Ez(2,i)+k3*(Ez(1,i)-Ez0(2,i));
    Ez0(end,i)=Ez(end-1,i)+k3*(Ez(end,i)-Ez0(end-1,i));
end 
   Ez=Ez0; 
% %АШ будет на 1/2 новее
% %аш использует уже посчитанное Е
% %ЦИКЛ и й к
    for i=2:nx-1
        
            Ez(i,1)=exp(-((i-nx/2+0.5)^2*dx^2)/w^2+1i*omega*tand(fi)*(i-nx/2+0.5)*dx)*exp(1i*(t)*dt*omega)*exp(-(t)*dt/tau);%без затухания
            %Hx(i,1)=(1/Z)*exp(-((i-n/2+0.5)^2*dx^2)/w^2+1i*omega*tand(fi)*(i-n/2+0.5)*dx)*exp(1i*(t-0.5)*dt*omega + 1i*omega*dy/2)*exp(-(t-0.5)*dt/tau);% затухания
        
    end

    Hx0(1:end, 1:tt) = Hx(1:end, 1:tt) + dH *(Ez(1:end-1, 1:tt) - Ez(1:end-1, 2:tt+1));
    Hy0(1:end, 1:tt) = Hy(1:end, 1:tt) + dH *(Ez(2:end, 1:tt) - Ez(1:end-1, 1:tt));
%     for i=1:nx-1
%     for j=1:tt %ny-1
%         % константы упрощения
%             Hx0(i,j)=Hx(i,j)+dH*(Ez(i,j)-Ez(i,j+1));
%             Hy0(i,j)=Hy(i,j)+dH*(Ez(i+1,j)-Ez(i,j));
%     end
% end

    Hx=Hx0;Hy=Hy0; 

  if (mod(t,99)==0)
      subplot(1,2,1)
      
    imagesc((abs(Ez)).^2)
    colormap('hot')
    %colorbar
    c=clock-clo;
    sec=c(6)+60*(c(6)<0);
    time=c(3)*86400+c(4)*3600+c(5)*60+c(6);
    speed=t/time*10;
    estimated=(tmax-t)/speed/6+0.99;
    title([' Step ' int2str(t) '. Time ' int2str(c(5)-(c(6)<0)) ':' int2str(sec) ' Speed ' int2str(speed) '/10  Estimated +' int2str(estimated) ' min']);
    subplot(1,2,2)
    spectrTemp=single(abs(fftshift(fft(fftshift(Ez)))));
    imagesc(spectrTemp(nx/2-150:nx/2+150,:))
    %clear spectrTemp;
    colormap('hot')
    %pause(0.003);
    drawnow
end
end
c=clock-clo;
sec=c(6)+60*(c(6)<0);
time=c(3)*86400+c(4)*3600+c(5)*60+c(6)




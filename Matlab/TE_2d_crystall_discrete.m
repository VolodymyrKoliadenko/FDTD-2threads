clc; clear all; %close all;
%clo=clock
dx=0.5*1e-7/1.4;
dy=dx;% �� �����, ����� ������������ ���
period=1.2e-6;
Q=1.5;
n=1;%ne ������������ � ������� ������� ���� =1
lambd=1064e-9;%1064
prodol=2*n*period^2/lambd/Q;
nx=1700;%600
ny=900;%1700
Imeem__Sloy= ny*dy/prodol*2
Ez=zeros(nx,ny)+1i*1e-100;Ez0=Ez;%� ������� �������, ��� ������ ����������. ���� ��������:)
Hy=zeros(nx-1,ny-1)+1i*1e-100;Hy0=Hy;
Hx=zeros(nx-1,ny-1)+1i*1e-100;Hx0=Hx;
e=ones(nx,ny);%epsilon
i=2;j=2;k=2;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5


omega=2*pi/lambd;%������� HERZ /3e8
 
dt=dx/2;%����� ���=�� ���� 0,2
tau=50000*dt;
%dz=1;
Z=376.7303;
s=dt/dx;  k3=(1-s)/(1+s);
w=15*1e-7;% 25
ii=1:ny;
%_____      E P S I L O N
 
m=(0.008*2*4);% ������� � 2���� ������ ������ ���� /��������� ��� ������
nachalo=3;
% �������������� ���������
for i=1:nx
    e(i,1:nachalo-1)=1;
    e(i,nachalo:end)=1+(m/2)*(1+sign(-0.1+cos(2*pi*(i-nx/2+0.5)*dx/period)*sin(2*pi*(ii(nachalo:end)-nachalo)*dy/prodol)));
end
e=(dt*Z/dx)./e;
dH=dt/dx/Z;
%e=(dt*Z/dx)./e;
fi=0;%asind(lambd/2/period);%gradus

clo=clock;
%i=2: ���� � � �
tmax=ny*2.2;
for t=1:tmax
    tt=min(t*s+10,ny-1);
    for i=2:nx-1
        
            Ez(i,1)=exp(-((i-nx/2+0.5)^2*dx^2)/w^2+1i*omega*tand(fi)*(i-nx/2+0.5)*dx)*exp(1i*(t-1)*dt*omega)*exp(-(t-1)*dt/tau);% ���������
            %Hx(i,1)=(1/Z)*exp(-((i-n/2+0.5)^2*dx^2)/w^2+1i*omega*tand(fi)*(i-n/2+0.5)*dx)*exp(1i*(t-0.5)*dt*omega + 1i*omega*dy/2)*exp(-(t-0.5)*dt/tau);% ���������
        
    end
for i=2:nx-1
    for j=2:tt %ny-1
        Ez0(i,j)=Ez(i,j)+e(i,j)*((Hx(i,j-1)-Hx(i,j)+Hy(i,j)-Hy(i-1,j)));
            
    end
end
for i=1:nx
    Ez0(i,end)= Ez(i,end-1)+k3*(Ez(i,end)-Ez0(i,end-1)); %GOOD ������ ����
end
for i=1:ny
    Ez0(1,i)= Ez(2,i)+k3*(Ez(1,i)-Ez0(2,i));
    Ez0(end,i)=Ez(end-1,i)+k3*(Ez(end,i)-Ez0(end-1,i));
end 
   Ez=Ez0; 
% %�� ����� �� 1/2 �����
% %�� ���������� ��� ����������� �
% %���� � � �
    for i=2:nx-1
        
            Ez(i,1)=exp(-((i-nx/2+0.5)^2*dx^2)/w^2+1i*omega*tand(fi)*(i-nx/2+0.5)*dx)*exp(1i*(t)*dt*omega)*exp(-(t)*dt/tau);%��� ���������
            %Hx(i,1)=(1/Z)*exp(-((i-n/2+0.5)^2*dx^2)/w^2+1i*omega*tand(fi)*(i-n/2+0.5)*dx)*exp(1i*(t-0.5)*dt*omega + 1i*omega*dy/2)*exp(-(t-0.5)*dt/tau);% ���������
        
    end
for i=1:nx-1
    for j=1:tt %ny-1
        % ��������� ���������
            Hx0(i,j)=Hx(i,j)+dH*(Ez(i,j)-Ez(i,j+1));
            Hy0(i,j)=Hy(i,j)+dH*(Ez(i+1,j)-Ez(i,j));
    end
end
Hx=Hx0;Hy=Hy0; 

  if (mod(t,10)==0)
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
    imagesc(spectrTemp(nx/2-50:nx/2+50,:))
    %clear spectrTemp;
    colormap('hot')
    drawnow
end
end




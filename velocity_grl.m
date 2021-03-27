
rw=1000; % density water
ri=910; % density ice
rs=2900; % density solid
rg=1.8; % density gas
kw=2.2e9; % bulk modulus water
ki=8.1e9; % bulk modulus ice
ks=80e9; % bulk modulus solid
kg=1e5; % bulk modulus gas 
gw=0; % shear modulus water
gi=3.7e9; % shear modulus ice
gs=40e9; % shear modulus solid
gg=0; %shear modulus gas

n=150;
f=linspace(0,30,n); % range of porosity
a=logspace(-1.5,-0.01,n); % range of aspect ratio
vsice=sqrt(gi/ri);
vssolid=sqrt(gs/rs);
ps=(3*ks-2*gs)/(6*ks+2*gs);

%f=0.1;
%a=0.1;
vvsi=zeros(n,n);
vvsw=zeros(n,n);
vvsg=zeros(n,n);
vvpi=zeros(n,n);
vvpw=zeros(n,n);
vvpg=zeros(n,n);

for i=1:n;
    for j=1:n;
        aa=a(i);
        ff=f(j)/100; % pore fraction from porosity

        [vsi,vpi]=wu(rs,ri,ks,ki,gs,gi,aa,ff);
        [vsw,vpw]=gassmann(rs,rw,ks,kw,gs,gw,aa,ff);
        [vsg,vpg]=wu(rs,rg,ks,kg,gs,gg,aa,ff);


        vvsg(i,j)=vsg;
        vvsi(i,j)=vsi;
        vvsw(i,j)=vsw;
        vvpg(i,j)=vpg;
        vvpi(i,j)=vpi;
        vvpw(i,j)=vpw;
        vpvsg(i,j)=vpg/vsg;
        vpvsi(i,j)=vpi/vsi;
        vpvsw(i,j)=vpw/vsw;
        
    end
end

[F,A]=meshgrid(f,a);


figure (1)
hold off
contourf(F,A,0.001*vvsg,[1.2,1.6,2.,2.4,2.8,3.2,3.6],'ShowText','on')
yticks([0.01,0.02,0.05,0.1,0.2,0.5])
set(gca,'YScale','log')
colormap(flipud(cool))
ax = gca;
ax.FontSize = 14; 
title('Vs gas (km/s)','FontSize', 18)
ylabel('aspect ratio','FontSize', 18);
xlabel('porosity (%)','FontSize', 18);
hold on
%contour(f,a,vvsg,[1700,2100],'LineWidth',4,'linecolor','k','ShowText','off')
%set(gca,'YScale','log')


%figure (2)
%contour(F,A,vvsw-vvsi,[0,-200,-400,-600,-800,-1000],'ShowText','on')
%set(gca,'YScale','log')
%title('Vs(water)-Vs(ice)')

figure (3)
contourf(F,A,vvsi/1000,'ShowText','on')
yticks([0.01,0.02,0.05,0.1,0.2,0.5])
set(gca,'YScale','log')
colormap(flipud(cool))
ax = gca;
ax.FontSize = 14; 
title('Vs ice (km/s)','FontSize', 18)
ylabel('aspect ratio','FontSize', 18);
xlabel('porosity (%)','FontSize', 18);

figure(4)
contourf(f,a,vvsw/1000,[1.2,1.6,2.0,2.4,2.8,3.2,3.6],'ShowText','on')
yticks([0.01,0.02,0.05,0.1,0.2,0.5])
set(gca,'YScale','log')
colormap(flipud(cool))
ax = gca;
ax.FontSize = 14; 
title('Vs (water)','FontSize', 18)
ylabel('aspect ratio','FontSize', 18);
xlabel('porosity (%)','FontSize', 18);

figure(5)
contourf(f,a,vvsg-vvsw,'ShowText','on')
yticks([0.01,0.02,0.05,0.1,0.2,0.5])
set(gca,'YScale','log')
colormap(flipud(cool))
ax = gca;
ax.FontSize = 14; 
title('Vs gas - Vs water (m/s)','FontSize', 20)
ylabel('aspect ratio','FontSize', 18);
xlabel('porosity (%)','FontSize', 18);

%contour(f,a,vsi-vsw,'ShowText','on')

figure(6)
contourf(f,a,vpvsg,[1.4,1.5,1.6,1.7,1.8,1.9,2.0],'ShowText','on')
set(gca,'YScale','log')
colormap(flipud(cool))
ax = gca;
ax.FontSize = 14; 
title('Vp/Vs gas','FontSize', 18)
ylabel('aspect ratio','FontSize', 18);
xlabel('porosity (%)','FontSize', 18);

figure(7)
contourf(f,a,vpvsi,[1.6,1.7,1.8,1.9,2.0],'ShowText','on')
set(gca,'YScale','log')
colormap(flipud(cool))
ax = gca;
ax.FontSize = 14; 
title('Vp/Vs ice','FontSize', 18)
ylabel('aspect ratio','FontSize', 18);
xlabel('porosity (%)','FontSize', 18);


figure(8)
contourf(f,a,vpvsw,[1.6,1.7,1.8,1.9,2.0],'ShowText','on')
set(gca,'YScale','log')
colormap(flipud(cool))
ax = gca;
ax.FontSize = 14; 
title('Vp/Vs water','FontSize', 18)
ylabel('aspect ratio','FontSize', 18);
xlabel('porosity (%)','FontSize', 18);

function [vs,vp]=wu(rm,ri,km,ki,gm,gi,alp,f)
% use formulas in Heap (2019)
% m for matrix, i for inclusion

% initial guess for moduli by linear interpolation
gv=gm*(1-f)+gi*f;
kv=km*(1-f)+ki*f;

phi=alp*(acos(alp)-alp*sqrt(1-alp*alp))/((1-alp*alp)^1.5);
g=alp*alp*(3*phi-2)/(1-alp*alp);

del=abs((gm-gv)/gm);
while(del>0.00000001 & gv>0)
    
 a=(gi/gv)-1;
 b=(1/3)*((ki/kv)-(gi/gv));
 r=(3*gv)/(3*kv+4*gv);
 
 f1=1+a*(1.5*(g+phi)-r*(1.5*g+2.5*phi-(4/3)));
 f2=1+a*(1+1.5*(g+phi)-0.5*r*(3*g+5*phi))+b*(3-4*r)+0.5*a*(a+3*b)*(3-4*r)*(g+phi-r*(g-phi+2*phi*phi));
 f3=1+0.5*a*(r*(2-phi)+((1+alp*alp)*g*(r-1)/(alp*alp)));
 f4=1+0.25*a*(3*phi+g-r*(g-phi));
 f5=a*(r*(g+phi-4/3)-g)+b*phi*(3-4*r);
 f6=1+a*(1+g-r*(g+phi))+b*(1-phi)*(3-4*r);
 f7=2+0.25*a*(9*phi+3*g-r*(5*phi+3*g))+b*phi*(3-4*r);
 f8=a*(1-2*r+0.5*g*(r-1)+0.5*phi*(5*r-3))+b*(1-phi)*(3-4*r);
 f9=a*(g*(r-1)-r*phi)+b*phi*(3-4*r);
 tiijj=3*f1/f2;
 dt=(2/f3)+(1/f4)+((f4*f5+f6*f7-f8*f9)/(f2*f4));
 
 % Wu (1966) self-consistent model
 gc=gm+(f*(gi-gm)*dt/5);
 kc=km+(f*(ki-km)*tiijj/3);
 
 % the number and (1-number) set how quickly to converge
 gv=0.9*gv+0.1*gc;
 kv=0.9*kv+0.1*kc;
 del=abs((gv-gc)/gv);
end

 rc=rm*(1-f)+f*ri;
 if gv>0
   vs=sqrt(gv/rc);
 else
   vs=NaN;
 end
 if gv>0 && kv > 0
   vp=sqrt((kv+(4/3)*gv)/rc);
 else
     vp=NaN;
 end
 
end



function [vs,vp]=gassmann(rm,ri,km,ki,gm,gi,alp,f)
% use Gassmann p 169 to saturate pores
% m for matrix i for inclusion
kdry=0;
gv=gm*(1-f)+gi*f;
kv=km*(1-f)+kdry*f;

phi=alp*(acos(alp)-alp*sqrt(1-alp*alp))/((1-alp*alp)^1.5);
g=alp*alp*(3*phi-2)/(1-alp*alp);

del=abs((gm-gv)/gm);
while(del>0.00000001 & gv>0)
    
 a=(gi/gv)-1;
 b=(1/3)*((kdry/kv)-(gi/gv));
 r=(3*gv)/(3*kv+4*gv);
 
 f1=1+a*(1.5*(g+phi)-r*(1.5*g+2.5*phi-(4/3)));
 f2=1+a*(1+1.5*(g+phi)-0.5*r*(3*g+5*phi))+b*(3-4*r)+0.5*a*(a+3*b)*(3-4*r)*(g+phi-r*(g-phi+2*phi*phi));
 f3=1+0.5*a*(r*(2-phi)+((1+alp*alp)*g*(r-1)/(alp*alp)));
 f4=1+0.25*a*(3*phi+g-r*(g-phi));
 f5=a*(r*(g+phi-4/3)-g)+b*phi*(3-4*r);
 f6=1+a*(1+g-r*(g+phi))+b*(1-phi)*(3-4*r);
 f7=2+0.25*a*(9*phi+3*g-r*(5*phi+3*g))+b*phi*(3-4*r);
 f8=a*(1-2*r+0.5*g*(r-1)+0.5*phi*(5*r-3))+b*(1-phi)*(3-4*r);
 f9=a*(g*(r-1)-r*phi)+b*phi*(3-4*r);
 tiijj=3*f1/f2;
 dt=(2/f3)+(1/f4)+((f4*f5+f6*f7-f8*f9)/(f2*f4));
 
 % Wu (1966) self-consistent model 
 gc=gm+(f*(gi-gm)*dt/5);
 kc=km+(f*(kdry-km)*tiijj/3);
 
 % the number and (1-number) set how quickly to converge
 gv=0.9*gv+0.1*gc;
 kv=0.9*kv+0.1*kc;
 del=abs((gv-gc)/gv);
end

 rc=rm*(1-f)+f*ri;
 % use gassmann, page 169 in Mavko et al. (1998)
 kv=kv+((1-(kv/km))^2)/((f/ki)+((1-f)/km)-(kv/(km*km)));
 
 if gv>0
   vs=sqrt(gv/rc);
 else
   vs=NaN;
 end
 if gv>0 && kv > 0
   vp=sqrt((kv+(4/3)*gv)/rc);
 else
     vp=NaN;
 end
 
end

function [vs,vp]=oconnell(rm,ri,km,gm,f)
% use formulas in Mavko
% m for matrix i for inclusion

pm=(3*km-2*gm)/(6*km+2*gm); % poisson ratio solid

psc=pm*(1-16*f/9);
ksc=km*(1-(16/9)*f*((1-psc*psc)/(1-2*psc)));
gsc=gm*(1-(32/45)*f*((1-psc)*(5-psc)/(2-psc)));

rc=rm*(1-f)+f*ri;
 
vs=sqrt(gsc/rc);
vp=sqrt((ksc+(4/3)*gsc)/rc);
 
 
end


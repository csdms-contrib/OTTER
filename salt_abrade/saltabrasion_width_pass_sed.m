close all
clear all

%Model based on Yanites, 2018 JGR manuscript: The dynamics of channel slope, width, and sediment in actively eroding bedrock river systems

%model set up
upinc=5
dx=2000

x=5000:dx:205000;
model_name=['saltabrade_']

%hack's law values calculated by regression profile 1, L vs. A (both in m)
Hc=1
He=1.8598
Ah=Hc.*(x.^(He));


% initialize fluvial variables
g=9.8;
rhow=1000;
rhos=2650;
kf=-10^-5;
thresh=0;
R=(rhos./rhow)-1;

lambda=0.2; %sediment porosity
n=0.04
tauc=0.0495; %from Wonga nd parker 2006



%uplift set-up
Uprate=0.5*(10^-3);
uplift=zeros(1,length(x))+Uprate;
% fold: uplift(10:15)=Uprate.*5


kQ=10^(-7)
eQ=1
Qw=kQ.*(Ah.^eQ);


kw=5
ew=0.5;
W=kw.*(Qw.^ew);
Wo=W;

%%sed supply
beta=0.3;
yr2sec=3.14.*(10^7);
onoff=1;
Qs=onoff.*rhos.*beta.*Ah.*uplift./yr2sec;
dQs=[Qs(1) diff(Qs)];
dA=[Ah(1) diff(Ah)];
%sediment characteriscs
%sternberg's law
a=0.02./(10^3);
%next line is for capture
%LA=[x(x<=Lcap.*1000) x(x>Lcap.*1000)-Lcap.*1000];

%the normal "distance downstream" is replaced by a crude hack's law for capture

D0=0.2
D=D0.*exp(-((Ah./Hc).^(1./He)).*a);
D(D<0.005)=0.005;


%topogrpahy set up
z=zeros(1,length(x))+10;
z=.007.*(max(x)-x)
sed_depth=zeros(1,length(x))+.1;
topo=z+sed_depth;
depth_threshold=0.2;
%time set up
time=20*(10^6);
dt=5;

tarray=1:dt:time;

%variable set up
Qs_down=0;
dz_s=zeros(1,length(x));
slope=zeros(1,length(x));
qs=zeros(1,length(x));
depth_threshold=D.*1;
dz_b=zeros(1,length(x));

dz_b_t=zeros(10000,length(dz_b));

eroded_sedsupply=uplift;
last_supply=eroded_sedsupply;
tau_b=zeros(1,length(x));
H=zeros(1,length(x));
F=zeros(1,length(x));
Qt=zeros(1,length(x));
dz_s=zeros(1,length(x));
F=zeros(1,length(x));
Eis=zeros(1,length(x));
dz_b_store=dz_b;
dz_s_store=dz_s;
%% saltation abrasion parameters
Rb=(rhos-rhow)./rhow;

Y=5.*(10^10);
sigmaT=14.*(10^6); %pascals

%sklar correctio to 2004 says 10^6,
kv=10^6;

%% calculate settling velocities
visc=.001./1000;
Dstar=(rhos-rhow).*g.*(D.^3)./(rhow.*(visc.^2));
wstar=zeros(1,length(Dstar));
lgW=-3.76715+(1.92944.*log(Dstar))-(0.09815.*((log(Dstar)).^2))-(.00575.*((log(Dstar)).^3))+(.00056.*((log(Dstar)).^4));
wstar(Dstar>0.05)=exp(lgW(Dstar>0.05));
wstar(Dstar<=0.05)=(1.71.*(10^-4)).*(D(Dstar<=0.05).^2);
wf=(wstar.*(rhos-rhow).*g.*visc./rhow).^(1./3);

prefac=0.08.*Rb.*g.*Y./(kv.*(sigmaT.^2));
frac_yr_transport=.1;
k=0
%% loop
onlyonce=0

for i=1:length(tarray)
    
    z=z+(uplift.*dt);
    z(end)=0;
    sed_depth(end)=0;
    topo=z+sed_depth;
    slope=[diff(topo)./diff(x) 0];
    %slope=[slope 0];
    %calc shear stress
    %     tau_b=rhow.*g.*((Qw./W).^0.6).*((-1.*slope).^0.7).*dt;
    %     tau_b(tau_b<0)=0;
    %crit shear stress to entrain
    if tarray(i)>5000000 && onlyonce==0
        uplift=uplift.*upinc;
        onlyonce=1;
    end
    
    if mod(i,1)==0 || i==1
        for j=1:length(x)
            if slope(j)<0
                
                Qs_n=(rhos.*beta.*dA(j).*(eroded_sedsupply(j)).*dt)+Qs_down;
                % Qs_n=(rhos.*beta.*dA(j).*(uplift(j)).*dt)+Qs_down;
                
                %% NEED NEW CHANNEL WIDTH FUNCTION FOR SALT-ABRASION MODEL
                
                [Wo(j),H(j),F(j),Qtout]=calc_width_saltabr(dz_b_store(j)./dt,Qs_n./(yr2sec.*dt),Qw(j),kf,-slope(j),W(j),rhos,D(j),tauc,wf(j),prefac,frac_yr_transport,dt,dz_s_store(j),H(j));
                
                
                if isreal(Qtout)==1
                    Qt(j)   = Qtout;
                else
                    Qt(j)=0;
                end
                
                dz_s(j)=(1-(1-lambda)).*((dt.*Qt(j).*yr2sec)-Qs_n)./(W(j).*rhos.*dx);
                
                Qs_pass=dt.*Qt(j).*yr2sec;%.*W(j);
                if dz_s(j)>sed_depth(j) && sed_depth(j)>=0
                    dz_s(j)=sed_depth(j);
                    Qs_pass=Qs_n;
                end
                Qs_down=(rhos.*dz_s(j).*W(j).*(dx))+(Qs_pass);
                Qs(j)=Qs_n./dt;
                Qs_n=0;
                
                if Qs_down<0
                    Qs_down=0;
                end
            end
        end
    end
    sed_depth=sed_depth-(dz_s);
    sed_depth(sed_depth<0)=0;
    
    W=Wo;
    tau_b=rhow.*g.*(-slope).*(W.*H./(W+(2.*H)));
    %% modify for saltation abrasion model
    
    tau_star=tau_b./((rhos-rhow).*g.*D);
    ustar=((W.*H./(W+(2.*H))).*g.*(-slope)).^0.5;
    %erosion in m/s
    Eis=prefac.*(Qs./W).*(((tau_star./tauc)-1).^(-0.52)).*F.*(((1-((ustar./wf).^(2))).^(3./2)));
    Eis(Eis<(10^-5) & Eis>0)=10^-5;
    %dz_b=-Eis.*(pi.*(10^7)).*dt;
    dz_b=-Eis.*dt;
    dz_b(dz_b>0)=0;
    dz_b(sed_depth>=depth_threshold)=0;
    if i<10000
        dz_b_t(i,:)=dz_b./dt;
        
    elseif i>=10000
        ff=mod(i,10000)+1;
        dz_b_t(ff,:)=dz_b./dt;
    end
    
    
    if i>10000
        
        %         dz_last=((1000.*dt).*uplift)-(z-last_z);
        %         eroded_sedsupply=dz_last./(1000.*dt);
        eroded_sedsupply= -1.*mean(dz_b_t);
        
        % eroded_sedsupply(eroded_sedsupply<(1*10^-5))=(1*10^-5);
        %%%CHECK THIS
        eroded_sedsupply(eroded_sedsupply<(0.1.*uplift))=0.1.*uplift(eroded_sedsupply<(0.1.*uplift));%last_supply(eroded_sedsupply<(1*10^-5));
        %         eroded_sedsupply(eroded_sedsupply>(3.*last_supply))=last_supply(eroded_sedsupply>(3.*last_supply)).*3;
        
        %eroded_sedsupply(eroded_sedsupply>uplift)=uplift(eroded_sedsupply>uplift);
        
        % eroded_sedsupply(eroded_sedsupply>(1.001.*last_supply))=last_supply(eroded_sedsupply>(1.001.*last_supply)).*1.001;
        % eroded_sedsupply(eroded_sedsupply<(0.999.*last_supply))=last_supply(eroded_sedsupply<(0.999.*last_supply)).*.999;
        
        
        
        last_supply=eroded_sedsupply;
        
        last_z=z;
    end
    
    z(sed_depth<depth_threshold & slope<0)=z(sed_depth<depth_threshold & slope<0)+dz_b(sed_depth<depth_threshold & slope<0);
    z(end)=0;
    sed_depth(end)=0;
    
    topo=z+sed_depth;
    dz_b_store=dz_b;
    dz_s_store=dz_s;
    if mod(i,50000./dt)==0 || i==1
        k=k+1;
        %        hold off
        %
        %        plot(x./1000,z+sed_depth,x./1000,z)
        %        hold on
        %        plot(x([min(MR1) max(MR1)])./1000,topo([min(MR1) max(MR1)]),'kx')
        %        plot(x([min(URG1)])./1000,topo([min(URG1)]),'gx')
        %        title([num2str(i.*dt)])
        %        drawnow
        %  W(end)=W(end-1);
        
        subplot1=subplot(5,2,1);
        plot(x./1000,topo,'b')
        hold on
        plot(x./1000,z,'k')
        
        ylabel('Elevation (m)')
        hold on
        
        title(['Time = ' num2str(i.*dt)])
        
        subplot2=subplot(5,2,3);
        plot(x./1000,W)
        ylabel('Width (m)')
        subplot3=subplot(5,2,5);
        plot(x(1:end)./1000,F(1:end),'r')
        % plot(x(1:end)./1000,sed_depth(1:end),'r')
        
        xlabel('Distance (km)')
        ylabel('F')
        subplot9=subplot(5,2,7);
        %        plot(x(1:end)./1000,Qt,'r')
        %        hold on
        %        plot(x(1:end)./1000,Qs(1:end)./yr2sec,'b')
        plot(x(1:end)./1000,eroded_sedsupply)
        xlabel('Distance (km)')
        ylabel('sed trans cap/flux')
        subplot4=subplot(5,2,9);
        plot(x(1:end)./1000,-1.*dz_b,'r')
        hold on
        plot(x(1:end)./1000,uplift,'k')
        xlabel('Distance (km)')
        ylabel('bedrock erosion')
        drawnow
        subplot5=subplot(5,2,2);
        scatter(Qw,W)
        xlabel('Water Discharge (m^3/s)')
        ylabel('Channel Width (m)')
        pQw=polyfit(log10(Qw(2:end-1)),log10(W(2:end-1)),1);
        title(['b= ' num2str(pQw(1))])
        
        
        subplot7=subplot(5,2,4);
        scatter(-slope,W)
        xlabel('Channel Slope')
        ylabel('Channel Width (m)')
        pSw=polyfit(log10(-slope(2:end-1)),log10(W(2:end-1)),1);
        title(['c= ' num2str(pSw(1))])
        
        subplot6=subplot(5,2,6);
        scatter(-dz_b./dt,W)
        xlabel('Erosion Rate')
        ylabel('Channel Width (m)')
        pEw=polyfit(log10(-dz_b(2:end-1)),log10(W(2:end-1)),1);
        title(['d= ' num2str(pEw(1))])
        
        subplot8=subplot(5,2,8);
        scatter(x./1000,D)
        ylabel('Grain Size')
        xlabel('Distance Downstream')
        subplot10=subplot(5,2,10);
        scatter(x./1000,W./H)
        ylabel('W/H ratio')
        xlabel('Distance Downstream')
        
        
        %                end
        drawnow
        fig = gcf;
        fig.PaperUnits = 'inches';
        fig.PaperPosition = [0 0 6 6];
        fig.PaperPositionMode = 'manual';
        print('-dpng',['./moviefiles/' model_name '_n_' num2str(k)])
        clear fig
        hold(subplot1,'off')
        hold(subplot2,'off')
        hold(subplot3,'off')
        hold(subplot4,'off')
        hold(subplot5,'off')
        hold(subplot6,'off')
        hold(subplot7,'off')
        hold(subplot8,'off')
        hold(subplot9,'off')
        hold(subplot10,'off')
        
        
        cla(subplot1)
        cla(subplot2)
        cla(subplot3)
        cla(subplot4)
        cla(subplot5)
        cla(subplot6)
        cla(subplot7)
        cla(subplot8)
        cla(subplot9)
        cla(subplot10)
        hold off
        
        clear subplot*
        save(['./output/' model_name '_n_' num2str(k)],'-regexp', '^(?!(tarray|dz_b_t)$).')
    end
    %reset variables
    Qs_down=0;
    dz_s(1:end)=0;
    dz_b(1:end)=0;
end

k=k+1;
%        hold off
%
%        plot(x./1000,z+sed_depth,x./1000,z)
%        hold on
%        plot(x([min(MR1) max(MR1)])./1000,topo([min(MR1) max(MR1)]),'kx')
%        plot(x([min(URG1)])./1000,topo([min(URG1)]),'gx')
%        title([num2str(i.*dt)])
%        drawnow
W(end)=W(end-1);

subplot1=subplot(4,1,1);
plot(x(1:end)./1000,topo(1:end),'b')
hold on
plot(x(1:end)./1000,z(1:end),'k')

ylabel('Elevation (m)')
hold on

title(['Time = ' num2str(i.*dt)])

subplot2=subplot(5,1,2);
plot(x(1:end)./1000,W(1:end))
ylabel('Width (m)')
subplot3=subplot(5,1,3);
plot(x(1:end)./1000,sed_depth(1:end),'r')
xlabel('Distance (km)')
ylabel('sed depth')
subplot4=subplot(5,1,4);
plot(x(1:end)./1000,Qt(1:end),'r')
xlabel('Distance (km)')
ylabel('trans cap')
subplot4=subplot(5,1,5);
plot(x(1:end)./1000,dz_b,'r')
xlabel('Distance (km)')
ylabel('bedrock erosion')
drawnow

hold(subplot1,'off')
hold(subplot2,'off')
hold(subplot3,'off')
hold(subplot4,'off')

save(['./output/' model_name '_n_' num2str(k)])
print('-dpng',['./moviefiles/' model_name '_n_' num2str(k)])




close all
clear all

%Model based on Yanites, 2018 JGR manuscript: The dynamics of channel slope, width, and sediment in actively eroding bedrock river systems

%model set up
dx=2000
upinc=5
x=5000:dx:205000;
model_name=['sedcover_']

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


kQ=10^(-7)
eQ=1
Qw=kQ.*(Ah.^eQ);


kw=2
ew=0.5;
W=kw.*(Qw.^0.5);
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
D0=0.2;
D=D0.*exp(-((Ah./Hc).^(1./He)).*a);
D(D<0.0005)=0.0005;
plot(x./1000,D)

%topogrpahy set up
z=zeros(1,length(x))+(max(x)-x).*.010;
sed_depth=zeros(1,length(x));
topo=z+sed_depth;
depth_threshold=0.2;
%time set up
time=20*(10^6);
dt=1;

tarray=1:dt:time;

%variable set up
Qs_down=0;
dz_s=zeros(1,length(x));
slope=zeros(1,length(x));
qs=zeros(1,length(x));
depth_threshold=D.*3;
dz_b=zeros(1,length(x));
tau_b=zeros(1,length(x));
H=zeros(1,length(x));
F=zeros(1,length(x));
Qt=zeros(1,length(x));
dz_s=zeros(1,length(x));
k=0
dz_b_store=dz_b;
dz_s_store=dz_s;
%what is the fraction of the year that active transport and incision occurs?
frac_yr_transport=0.1;
dooonce=0
%% loop
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
    
    if tarray(i)>5000000 && dooonce==0
        doonce=1;
        uplift(:)=Uprate.*upinc;
    end
    %%MODIDFY TO PASS SEDIMENT IN DOWNSTREAM DIRECTION
    for j=1:length(x)
        if slope(j)<0 %% Do we have a river?
            
            
            Qs_n=0;
            
            %call calc_width for optimization of channel
            %geometry
            [Wo(j),H(j),DW]=calc_width(dz_b_store(j)./dt,Qs_n./(yr2sec.*dt),Qw(j),kf,-slope(j),W(j),rhos,D(j),tauc,frac_yr_transport,dt,dz_s_store(j),H(j));
            
        end
        
    end
    
    W=Wo;
    tau_b=rhow.*g.*(-slope).*(W.*H./(W+(2.*H)));
    
    dz_b=kf.*tau_b.*dt;
    dz_b(dz_b>0)=0;
    
    z=z+dz_b;
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
        W(end)=W(end-1);
        
        subplot1=subplot(5,2,1);
        plot(x(2:end)./1000,topo(2:end),'b')
        hold on
        plot(x(2:end)./1000,z(2:end),'k')
        
        ylabel('Elevation (m)')
        hold on
        
        title(['Time = ' num2str(i.*dt)])
        
        subplot2=subplot(5,2,3);
        plot(x(2:end)./1000,W(2:end))
        ylabel('Width (m)')
        subplot3=subplot(5,2,5);
        plot(x(1:end)./1000,F(1:end),'r')
        % plot(x(1:end)./1000,sed_depth(1:end),'r')
        
        xlabel('Distance (km)')
        ylabel('F')
        subplot9=subplot(5,2,7);
        plot(x(1:end)./1000,Qt,'r')
        hold on
        plot(x(1:end)./1000,Qs(1:end)./yr2sec,'b')
        xlabel('Distance (km)')
        ylabel('sed trans cap/flux')
        subplot4=subplot(5,2,9);
        plot(x(1:end)./1000,dz_b,'r')
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
        
        
        save(['./output/' model_name '_n_' num2str(k)])
        
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
    end
    %reset variables
    Qs_down=0;
    dz_s(1:end)=0;
    dz_b(1:end)=0;
end

beep




beep
%
%        subplot1=subplot(3,1,1);
%        plot(x(2:end),z(2:end),'b')
%        ylabel('Elevation (m)')
%        hold on
%        subplot2=subplot(3,1,2);
%        plot(x(2:end),W(2:end))
%        ylabel('Width (m)')
%        subplot3=subplot(3,1,3);
%        plot(x(2:end),tau_b(2:end),'r')
%        xlabel('Distance (km)')
%        ylabel('shear stress')
%


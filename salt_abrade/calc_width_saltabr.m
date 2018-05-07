function [Wo,H,F,Qto] = calc_width_saltabr(EE,QSS,Q,kv,S,Wi,rhos,D,tausc,wf,prefac,frac_yr_transport,dt,DS,Hin)
% calc_width_saltabr is a function that updates a given channel width based on the
% user input parameters and the optimizatoin scheme of Yanites and Tucker,
% 2010/turowski et al 2007 (but with Sklar cover model). Essentially, it calculates the current shear stress (no wide channel assumption) and
% the shear stress if channel width is slighlty larger.  Using these
% values, erosion potential is calculated for all thre scenarios (no change
% in width, slightly wider, slightly narrower) and the one that generates
% the maximum erosion rate dictatest the direction of channel change

%MODIFY FOR SKLAR MODEL
%  EE   erosion rate (not needed)
%  QSS  Sediment supply
%  Q    water discharge
%  kv   bedrock erodibility (not needed)
%  S    current slope
%  Wi   current width
%  rhos  sediment density
%  D    sediment size
% tausc critical shear stress
%
% outputs a new channel width,Wo, as a step increase or decrease (user define
% step percentage).
F=1;
rho=1000;
g=9.8;
n=.04;%manning's roughness
W=Wi;%logspace(log10(Wi)-1,log10(Wi)+1,1000);%3:.5:50;50;%W=(S.^(-3./16)).*(Q.^(3./8)).*(n.^3./8);  Channel Width array
SA=S;
%incase width is optimized.
Wo=W;
Wbig=W.*1.1;
Wsmall=W.*0.9;
Win=[Wi Wbig Wsmall];

S=SA;
%options=optimset('Display','off');
if W>5.*Hin
    Hout = [((n.*Q./W).^(3./5)).*(S.^(-3./10)) ((n.*Q./Wbig).^(3./5)).*(S.^(-3./10)) ((n.*Q./Wsmall).^(3./5)).*(S.^(-3./10))];
elseif W<=5.*Hin
    
    Htest= [((n.*Q./W).^(3./5)).*(S.^(-3./10)) ((n.*Q./Wbig).^(3./5)).*(S.^(-3./10)) ((n.*Q./Wsmall).^(3./5)).*(S.^(-3./10))];
    Hout=fsolve(@(Hw)((1./n).*Win.*Hw.*((Win.*Hw./((2.*Hw)+Win)).^(2./3)).*(S.^(1./2)))-Q,Htest,optimset('Display','off'));%,optimset('Display','off','Algorithm', {'levenberg-marquardt',10^-30}));%,'MaxIter',5) ); %%%find proper flow depth
    % [Hout ssq cnt]=LMFsolve(@(Hw)((1./n).*Win.*Hw.*((Win.*Hw./((2.*Hw)+Win)).^(2./3)).*(S.^(1./2)))-Q,Htest,'Display','0');%,'MaxIter',5) ); %%%find proper flow depth
    
end

%Hout=lsqnonlin(@(Hw)((1./n).*Win.*Hw.*((Win.*Hw./((2.*Hw)+Win)).^(2./3)).*
%(S.^(1./2)))-Q,Htest,Htest.*0.1,Htest.*1.4);%[0
%0],5000);,'NonlEqnAlgorithm', 'gn' ,'LargeScale','off'
H=Hout(1);
Hbig=Hout(2);
Hsmall=Hout(3);
%
% if isnan(H)==1
%     H=0;
% end
%Width=W;
% Wbig=W.*1.1;
% Htest= ((n.*Q./Wbig).^(3./5)).*(S.^(-3./10));
%
% Hbig=fsolve(@(Hw)((1./n).*Wbig.*Hw.*((Wbig.*Hw./((2.*Hw)+Wbig)).^(2./3)).*(S.^(1./2)))-Q,Htest,optimset('Display','off'));  %%%find proper flow depth
% % if isnan(Hbig)==1
%     Hbig=0;
% end
taunow=rho.*g.*S.*(W.*H./(W+(2.*H)));
taubigW=rho.*g.*S.*(Wbig.*Hbig./(Wbig+(2.*Hbig)));
tausmallW=rho.*g.*S.*(Wsmall.*Hsmall./(Wsmall+(2.*Hsmall)));

%     tausmW=rho.*g.*S.*(0.9.*W.*H./(0.9.*W+(2.*H)))


tausF=taunow./((rhos-rho).*g.*D); %%non-d shear
tausFbig=taubigW./((rhos-rho).*g.*D); %%non-d shear
tausFsmall=tausmallW./((rhos-rho).*g.*D); %%non-d shear

ustarF=((W.*H./(W+(2.*H))).*g.*(S)).^0.5;
ustarFbig=((Wbig.*Hbig./(Wbig+(2.*Hbig))).*g.*(S)).^0.5;
ustarFsmall=((Wsmall.*Hsmall./(Wsmall+(2.*Hsmall))).*g.*(S)).^0.5;






qssf=3.97.*((tausF-tausc).^(3./2)); %% Wong and Parker reanalaysis of MPM
qssfbig=3.97.*((tausFbig-tausc).^(3./2)); %% Wong and Parker reanalaysis of MPM
qssfsmall=3.97.*((tausFsmall-tausc).^(3./2)); %% Wong and Parker reanalaysis of MPM

qsf=sqrt((rhos-rho).*g.*D./rho).*D.*qssf; % sed trans per unit width, with wall
qsfbig=sqrt((rhos-rho).*g.*D./rho).*D.*qssfbig; % sed trans per unit width, with wall
qsfsmall=sqrt((rhos-rho).*g.*D./rho).*D.*qssfsmall; % sed trans per unit width, with wall

if qsf<0
    qsf=0;
end
if qsfbig<0
    qsfbig=0;
end
if qsfsmall<0
    qsfsmall=0;
end


Qtf=frac_yr_transport.*rhos.*qsf.*W; %% mass tranposrt capacity kg/s- with wall effects
Qto=Qtf;
Qtfbig=frac_yr_transport.*rhos.*qsfbig.*Wbig; %% mass tranposrt capacity kg/s- with wall effects
Qtfsmall=frac_yr_transport.*rhos.*qsfsmall.*Wsmall; %% mass tranposrt capacity kg/s- with wall effects

% if Qtf<=QAA
%     QSS=Qtf;
% end



F=(1-(QSS./Qtf));
Fbig=(1-(QSS./Qtfbig));
Fsmall=(1-(QSS./Qtfsmall));

%tau_star=tau_b./((rhos-rhow).*g.*D);

Eis=prefac.*((QSS./rhos)./W).*(((tausF./tausc)-1).^(-0.52)).*F.*(((1-((ustarF./wf).^(2))).^(3./2)));
Eisbig=prefac.*((QSS./rhos)./Wbig).*(((tausFbig./tausc)-1).^(-0.52)).*Fbig.*(((1-((ustarFbig./wf).^(2))).^(3./2)));
Eissmall=prefac.*((QSS./rhos)./Wsmall).*(((tausFsmall./tausc)-1).^(-0.52)).*Fsmall.*(((1-((ustarFsmall./wf).^(2))).^(3./2)));




% if Qtf<=QAA
%     QSS=Qtf;
% end
taubQSQT=taunow.*(1-(QSS./Qtf));
taubQSQTbig=taubigW.*(1-(QSS./Qtfbig));
taubQSQTsmall=tausmallW.*(1-(QSS./Qtfsmall));


Narrowing_rate_rock=W./(H./(EE))./50;
Narrowing_rate_sed=W./(H./(DS))./50;

% if Narrowing_rate_rock<Narrowing_rate_sed && Narrowing_rate_rock~=0
%
% elseif Narrowing_rate_rock>=Narrowing_rate_sed && Narrowing_rate_sed~=0
%
%
% elseif
Narrowing_rate=Narrowing_rate_rock;

if Narrowing_rate>0
    Narrowing_rate=Narrowing_rate.*(-1);
end

% limit narrowing to 10% of current width
if (Narrowing_rate)<(-.1.*Wi)
    Narrowing_rate=-(0.1.*Wi);
end

Kw=.0002;
Widening=Kw.*taunow;
%oversimplification need to fix when dynamically evolving sediment
% if QSS>Qtf
%     F=0;
% end

if QSS<=0
    
    if tausmallW>taubigW & tausmallW>taunow
        Wo=Wi+(dt.*Narrowing_rate);
        F=1;
    elseif taubigW>taunow & taubigW>tausmallW
        Wo =(dt.*Widening)+Wi;
        F=1;
    elseif taunow>=taubigW & taunow>=tausmallW
        Wo=Wi;
        F=1;
    end
    
elseif Qtf>=QSS & QSS>0
    if Eissmall>Eis & Eissmall>Eisbig
        Wo=Wi+(dt.*Narrowing_rate);
        
        %cover effect
        F=(1-(QSS./Qtfsmall));
        Qto=Qtfsmall;
        
    elseif Eisbig>Eis & Eisbig>Eissmall
        Wo =(dt.*Widening)+Wi;
        %cover effect
        Qto=Qtfbig;
        F=(1-(QSS./Qtfbig));
        
    elseif Eis>=Eissmall & Eis>=Eisbig
        Wo=Wi;
        Qto=Qtf;
        F=(1-(QSS./Qtf));
        
    end
    
    
    
    
elseif Qtf<QSS & QSS>0
    
    F=0;
    if Qtfsmall>Qtf & Qtfsmall>Qtfbig
        Wo=Wi+(dt.*Narrowing_rate);
        Qto=Qtfsmall;
    elseif Qtfbig>Qtf & Qtfbig>Qtfsmall
        Wo =(dt.*Widening)+Wi;
        Qto=Qtfbig;
    elseif Qtf>=Qtfsmall & Qtf>=Qtfbig
        Wo=Wi;
        Qto=Qtf;
    end
    
    
    
end

if isnan(F)==1
    F=0;
end

end
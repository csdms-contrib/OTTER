
function [Wo,H,DW] = calc_width(EE,QSS,Q,kv,S,Wi,rhos,D,tausc,frac_yr_transport,dt,DS,Hin)
% calc_width is a function that updates a given channel width based on the
% user input parameters and the optimizatoin scheme of Yanites and Tucker,
% 2009. Essentially, it calculates the current shear stress (no wide channel assumption) and
% the shear stress if channel width is slighlty larger.  If shear stress increases
% with an increase in width, channel width is incrementally increased.  If sheare
% stress decreases with an increase in width, channel width is decreased.  So the point is to
% always increase shear stress, driving the system towards optimatilly. Inputs are
%  EE   erosion rate (not needed)
%  QSS  Sediment supply
%  Q    water discharge
%  kv   bedrock erdability (not needed)
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
Wbig=W.*1.01;
Wsmall=W.*0.99;

Win=[Wi Wbig Wsmall];
S=SA;
%options=optimset('Display','off');

%set conditional for wide channel assumption
if W>5.*Hin
    Hout = [((n.*Q./W).^(3./5)).*(S.^(-3./10)) ((n.*Q./Wbig).^(3./5)).*(S.^(-3./10)) ((n.*Q./Wsmall).^(3./5)).*(S.^(-3./10))];
elseif W<=5.*Hin
    Htest= [((n.*Q./W).^(3./5)).*(S.^(-3./10)) ((n.*Q./Wbig).^(3./5)).*(S.^(-3./10)) ((n.*Q./Wsmall).^(3./5)).*(S.^(-3./10))];
    Hout=fsolve(@(Hw)((1./n).*Win.*Hw.*((Win.*Hw./((2.*Hw)+Win)).^(2./3)).*(S.^(1./2)))-Q,Htest,optimset('Display','off','Algorithm', {'levenberg-marquardt',10^-30}));%,'MaxIter',5) ); %%%find proper flow depth
    
    % Hout=fminunc(@(Hw)((1./n).*Win.*Hw.*((Win.*Hw./((2.*Hw)+Win)).^(2./3)).*(S.^(1./2)))-Q,Htest);%,'MaxIter',5) ); %%%find proper flow depth
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




%oversimplification need to fix when dynamically evolving sediment
% if QSS>Qtf
%     F=0;
% end

Narrowing_rate_rock=W./(H./(dt.*EE))./50;

Narrowing_rate_sed=W./(H./(dt.*DS))./50;


% if Narrowing_rate_rock<Narrowing_rate_sed && Narrowing_rate_rock~=0
%
% elseif Narrowing_rate_rock>=Narrowing_rate_sed && Narrowing_rate_sed~=0
%
%
% elseif
Narrowing_rate=min([Narrowing_rate_rock Narrowing_rate_sed]);

if Narrowing_rate>0
    Narrowing_rate=Narrowing_rate.*(-1);
end

% limit narrowing to 10% of current width
if (Narrowing_rate)<(-0.05.*Wi)
    Narrowing_rate=-(0.05.*Wi);
end

Kw=.0005;
Widening=Kw.*taunow;




if tausmallW>taubigW & tausmallW>taunow
    Wo=Wi+(dt.*Narrowing_rate);
elseif taubigW>taunow & taubigW>tausmallW
    Wo =(dt.*Widening)+Wi;
elseif taunow>=taubigW & taunow>=tausmallW
    Wo=Wi;
end





DW=Wi-Wo;
end

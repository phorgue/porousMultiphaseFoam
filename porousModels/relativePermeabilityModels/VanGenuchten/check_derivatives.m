% To check Van genuchten derivatives

%% symbolic model
syms S;

syms Smin;
syms Smax;

syms Se(S); %normalized saturation

Se(S) = (S-Smin) / (Smax - Smin);

syms m;%exposant

syms kra(S);
syms krb(S);



kra(S) = ((1-Se)^0.5) * (1-Se^(1/m))^2*m;
krb(S) = (Se^0.5) * (1-(1-(Se^(1/m)))^m)^2;

figure(1)
title('relative permeabilities VanGenuchten s Model')
m = 1;
Smin = 0.001;
Smax = .999;
ezplot(eval(kra),[0,1]);
hold on 
ezplot(eval(krb),[0,1]);
hold off


%% derivatives plotting
figure(2)
m = 1
Smin = 0.001;
Smax = .999;


syms dkra(S);
syms dkrb(S);

dkra(S) = diff(kra);
dkrb(S) = diff(krb);
% 
% ezplot(eval(dkra),[0,1]);
% hold on 
ezplot(eval(dkrb),[0,1])
hold off

%% calculated funtions to identify

figure(3)
syms dkra_c(S);
syms dkrb_c(S);

dkra_c(S) = - (1-(Se^(1/m)))^(2*m-1) * (-5*(Se^(1/m+1))+4*(Se^(1/m))+Se);
dkra_c(S) = dkra_c * 1/(2*((1-Se)^0.5)*Se);
dkra_c(S) = dkra_c * 1/(Smax - Smin);

dkrb_c(S) = 0.5 * (1-(1-(Se^(1/m))^m));
dkrb_c(S) = dkrb_c *( 4 * (Se^(1/m-1/2)) * ((1-(Se^(1/m)))^(m-1))) - (((1-(Se^(1/m)))^m) -1) / (Se^0.5);
dkrb_c(S) = dkrb_c * 1/(Smax - Smin)


% ezplot(eval(dkra_c),[0,1]);
% hold on 
ezplot(eval(dkrb),[0,1])
hold off
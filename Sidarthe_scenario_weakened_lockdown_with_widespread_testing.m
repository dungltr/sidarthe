%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MATLAB Code for epidemic simulations with the SIDARTHE model in the work
%
% Modelling the COVID-19 epidemic and implementation of population-wide interventions in Italy
% by Giulia Giordano, Franco Blanchini, Raffaele Bruno, Patrizio Colaneri, Alessandro Di Filippo, Angela Di Matteo, Marta Colaneri
% 
% Giulia Giordano, April 5, 2020
% Contact: giulia.giordano@unitn.it
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Italian population
popolazione=60e6;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Simulation horizon: CAN BE MODIFIED AT ONE'S WILL
Orizzonte = 350;

% Plot yes/no: SET TO 1 IF PDF FIGURES MUST BE GENERATED, 0 OTHERWISE
plotPDF = 0;

% Time-step for Euler discretisation of the continuous-time system
step=0.01;

% Transmission rate due to contacts with UNDETECTED asymptomatic infected
alfa=0.57;
% Transmission rate due to contacts with DETECTED asymptomatic infected
beta=0.0114;
% Transmission rate due to contacts with UNDETECTED symptomatic infected
gamma=0.456;
% Transmission rate due to contacts with DETECTED symptomatic infected
delta=0.0114;

% Detection rate for ASYMPTOMATIC
epsilon=0.171;
% Detection rate for SYMPTOMATIC
theta=0.3705;

% Worsening rate: UNDETECTED asymptomatic infected becomes symptomatic
zeta=0.1254;
% Worsening rate: DETECTED asymptomatic infected becomes symptomatic
eta=0.1254;

% Worsening rate: UNDETECTED symptomatic infected develop life-threatening
% symptoms
mu=0.0171;
% Worsening rate: DETECTED symptomatic infected develop life-threatening
% symptoms
nu=0.0274;

% Mortality rate for infected with life-threatening symptoms
tau=0.01;

% Recovery rate for undetected asymptomatic infected
lambda=0.0342;
% Recovery rate for detected asymptomatic infected
rho=0.0342;
% Recovery rate for undetected symptomatic infected
kappa=0.0171;
% Recovery rate for detected symptomatic infected
xi=0.0171;
% Recovery rate for life-threatened symptomatic infected
sigma=0.0171;

% Rate of loss of immunity
ki=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINITIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parameters
[r1,r2,r3,r4,r5] = calParamater(epsilon,zeta,lambda,eta,rho,theta,mu,kappa,nu,xi,sigma,tau)

% Initial R0
R0_iniziale = calculate1(alfa,r1,beta,epsilon,r2,gamma,zeta,r3,delta,eta,r4,theta)


% Time horizon
t=1:step:Orizzonte;

% Vectors for time evolution of variables
S=zeros(1,length(t));
I=zeros(1,length(t));
D=zeros(1,length(t));
A=zeros(1,length(t));
R=zeros(1,length(t));
T=zeros(1,length(t));
H=zeros(1,length(t));
H_diagnosticati=zeros(1,length(t)); % DIAGNOSED recovered only!
E=zeros(1,length(t));

% Vectors for time evolution of actual/perceived Case Fatality Rate
M=zeros(1,length(t));
P=zeros(1,length(t));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIAL CONDITIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

I(1)=200/popolazione;
D(1)=20/popolazione;
A(1)=1/popolazione;
R(1)=2/popolazione;
T(1)=0.00;
H(1)=0.00;
E(1)=0.00;
S(1)=1-I(1)-D(1)-A(1)-R(1)-T(1)-H(1)-E(1);

H_diagnosticati(1) = 0.00; % DIAGNOSED recovered only
Infetti_reali(1)=I(1)+D(1)+A(1)+R(1)+T(1); % Actual currently infected

M(1)=0;
P(1)=0;

% Whole state vector
x=[S(1);I(1);D(1);A(1);R(1);T(1);H(1);E(1);H_diagnosticati(1);Infetti_reali(1)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SIMULATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% "Control" binary variables to compute the new R0 every time a policy has
% changed the parameters
plottato = 0;
plottato1 = 0;
plottato_bis = 0;
plottato_tris = 0;
plottato_quat = 0;
plottato_fin = 0;

for i=2:length(t)
    
    if (i>4/step) % Basic social distancing (awareness, schools closed)
        [alfa, gamma, beta, delta] = Day4()
        if plottato == 0 % Compute the new R0
            [r1,r2,r3,r4,r5] = calParamater(epsilon,zeta,lambda,eta,rho,theta,mu,kappa,nu,xi,sigma,tau)
            R0_primemisure = calculate1(alfa,r1,beta,epsilon,r2,gamma,zeta,r3,delta,eta,r4,theta)
            plottato = 1;
        end
    end
    
    if (i>12/step)
        % Screening limited to / focused on symptomatic subjects
        epsilon = Day12()
        if plottato1 == 0
            [r1,r2,r3,r4,r5] = calParamater(epsilon,zeta,lambda,eta,rho,theta,mu,kappa,nu,xi,sigma,tau)
            R0_primemisureeps = calculate1(alfa,r1,beta,epsilon,r2,gamma,zeta,r3,delta,eta,r4,theta)
            plottato1 = 1;
        end
    end
    
    if (i>22/step) % Social distancing: lockdown, mild effect  
        [alfa, beta, gamma, delta, mu, nu, zeta, eta, lambda, rho, kappa, xi, sigma] = Day22()     
        if plottato_bis == 0 % Compute the new R0
            [r1,r2,r3,r4,r5] = calParamater(epsilon,zeta,lambda,eta,rho,theta,mu,kappa,nu,xi,sigma,tau)
            R0_secondemisure = calculate2(alfa,r1,beta,epsilon,r2,gamma,zeta,r3,delta,eta,r4,theta)
            plottato_bis = 1;
        end
    end
    
    if (i>28/step) % Social distancing: lockdown, strong effect
        [alfa, gamma] = Day28()
        if plottato_tris == 0 % Compute the new R0
            [r1,r2,r3,r4,r5] = calParamater(epsilon,zeta,lambda,eta,rho,theta,mu,kappa,nu,xi,sigma,tau)
            R0_terzemisure = calculate2(alfa,r1,beta,epsilon,r2,gamma,zeta,r3,delta,eta,r4,theta)
            plottato_tris = 1;
        end
    end
    
    
    if (i>38/step) % Broader diagnosis campaign
       [epsilon, rho, kappa, xi, sigma, zeta, eta] = Day38()
        if plottato_quat == 0 % Compute the new R0
            [r1,r2,r3,r4,r5] = calParamater(epsilon,zeta,lambda,eta,rho,theta,mu,kappa,nu,xi,sigma,tau)
            R0_quartemisure = calculate2(alfa,r1,beta,epsilon,r2,gamma,zeta,r3,delta,eta,r4,theta)
            plottato_quat = 1;
        end
    end
    
    % OUTLOOK: to EXPLORE POSSIBLE SCENARIOS,
    % change the parameters (after 50 days)
    if (i>50/step)
        [alfa, beta, gamma, delta, epsilon, theta, zeta, eta, mu, nu, tau, lambda, rho, kappa, xi, sigma, ki] = Day50WeekendLockdownWidespread()
        if plottato_fin == 0 % Compute the new R0     
            [r1,r2,r3,r4,r5] = calParamater(epsilon,zeta,lambda,eta,rho,theta,mu,kappa,nu,xi,sigma,tau)
            R0_final = calculate2(alfa,r1,beta,epsilon,r2,gamma,zeta,r3,delta,eta,r4,theta)
            plottato_fin = 1;
        end        
    end
    
    % Compute the system evolution
    
    B=[-alfa*x(2)-beta*x(3)-gamma*x(4)-delta*x(5) 0 0 0 0 0 ki 0 0 0;
        alfa*x(2)+beta*x(3)+gamma*x(4)+delta*x(5) -(epsilon+zeta+lambda) 0 0 0 0 0 0 0 0;
        0 epsilon  -(eta+rho) 0 0 0 0 0 0 0;
        0 zeta 0 -(theta+mu+kappa) 0 0 0 0 0 0;
        0 0 eta theta -(nu+xi) 0 0 0 0 0;
        0 0 0 mu nu  -(sigma+tau) 0 0 0 0;
        0 lambda rho kappa xi sigma -ki 0 0 0;
        0 0 0 0 0 tau 0 0 0 0;
        0 0 rho 0 xi sigma 0 0 0 0;
        alfa*x(2)+beta*x(3)+gamma*x(4)+delta*x(5) 0 0 0 0 0 0 0 0 0];
    x=x+B*x*step;
    
    % Update variables
    
    S(i)=x(1);
    I(i)=x(2);
    D(i)=x(3);
    A(i)=x(4);
    R(i)=x(5);
    T(i)=x(6);
    H(i)=x(7);
    E(i)=x(8);
    
    H_diagnosticati(i)=x(9);
    Infetti_reali(i)=x(10);
    
    % Update Case Fatality Rate
    
    M(i)=E(i)/(S(1)-S(i));
    P(i)=E(i)/((epsilon*r3+(theta+mu)*zeta)*(I(1)+S(1)-I(i)-S(i))/(r1*r3)+(theta+mu)*(A(1)-A(i))/r3);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FINAL VALUES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Variables
Sbar=S(length(t));
Ibar=I(length(t));
Dbar=D(length(t));
Abar=A(length(t));
Rbar=R(length(t));
Tbar=T(length(t));
Hbar=H(length(t));
Ebar=E(length(t));

% Case fatality rate
Mbar=M(length(t));
Pbar=P(length(t));

% Case fatality rate from formulas
Mbar1=Ebar/(S(1)-Sbar);
Pbar1=Ebar/((epsilon*r3+(theta+mu)*zeta)*(I(1)+S(1)-Sbar-Ibar)/(r1*r3)+(theta+mu)*(A(1)-Abar)/r3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
plot(t,Infetti_reali,'b',t,I+D+A+R+T,'r',t,H,'g',t,E,'k')
hold on
plot(t,D+R+T+E+H_diagnosticati,'--b',t,D+R+T,'--r',t,H_diagnosticati,'--g')
xlim([t(1) t(end)])
ylim([0 0.015])
%title('Actual vs. Diagnosed Epidemic Evolution')
xlabel('Time (days)')
ylabel('Cases (fraction of the population)')
legend({'Cumulative Infected','Current Total Infected', 'Recovered', 'Deaths','Diagnosed Cumulative Infected','Diagnosed Current Total Infected', 'Diagnosed Recovered'},'Location','east')
grid

if plotPDF==1
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperPosition', [0 0 24 16]);
    set(gcf, 'PaperSize', [24 16]); % dimension on x axis and y axis resp.
    print(gcf,'-dpdf', ['PanoramicaEpidemiaRealevsPercepita.pdf'])
end
%

figure
plot(t,I,'b',t,D,'c',t,A,'g',t,R,'m',t,T,'r')
xlim([t(1) t(end)])
ylim([0 1.1e-3])
%title('Infected, different stages, Diagnosed vs. Non Diagnosed')
xlabel('Time (days)')
ylabel('Cases (fraction of the population)')
legend({'Infected ND AS', 'Infected D AS', 'Infected ND S', 'Infected D S', 'Infected D IC'},'Location','northeast')
grid

if plotPDF==1
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperPosition', [0 0 24 16]);
    set(gcf, 'PaperSize', [24 16]); % dimension on x axis and y axis resp.
    print(gcf,'-dpdf', ['SuddivisioneInfetti.pdf'])
end
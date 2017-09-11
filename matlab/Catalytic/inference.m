%program to estimate dynamic transmission, age structured model from
%serological data
addpath(genpath('C:/Users/kamg2/Documents/takeda-norovirus'))
clear chain
Lmax=76;

    %load('serologicalPelosi')                        %load Italian serological data
load('serologicalSon')                           %load Korean serological data
    %N=58652875;                                 %Population size of Italy
N=50617045;                                      %Population size of Republic of Korea
age=serologicalSon(:,1);                                    %ages from serological data
positive=serologicalSon(:,3);                               %numbers positive in serological data
sample_size=serologicalSon(:,2);                            %sample size from serological data

%age specific mortality rates estimated for Italy and Korea using a
%Gompertz model
a=1:120;
% Italy(1)=0.000071797477560;                    
% Italy(2)=0.078331922874846;
% mu_whole=Italy(1)*exp(Italy(2)*a);
KO(1)=0.000048680478532;
KO(2)=0.081116006138959;
mu_whole=KO(1)*exp(KO(2)*a);

syms q nu delta

%age specific contact patterns
%load('ContactMatrixItaly')
% C=CR_ave_Italy; %
load('KoreanContact')
C=ones(Lmax);   %KoreanContact; 


%fixed parameters
gamma=1/2;          %recover rate
alpha=1;            %rate individuals move from latent to infectious
                    %age dependant death rate
mu=mu_whole(1:75)/365;
mu(76)=[-log(mean(exp(-cumsum(mu_whole(76:end)))))-sum(mu_whole(1:75)) ]/365;

%Build model equations
[Equations] = MakeMSEIRSequations(q,C, mu, nu, gamma , delta, alpha );  
%symbolic variables
M=sym('M',[Lmax,1]);S=sym('S',[Lmax,1]);E=sym('E',[Lmax,1]);I=sym('I',[Lmax,1]);R=sym('R',[Lmax,1]);
%ode variable
Y=sym('Y',[5*Lmax,1]);
%time
syms t
%sub in ode variable
EQNS=subs(Equations,[M.',S.',E.',I.',R.'],Y.');

%ode file
funcName='ODEfile' ;
%use matlabFunction to prepare optimised ode function
matlabFunction(EQNS, 'file', funcName,'vars', {t, Y,q,delta,nu});
disp('ready to start MCMC')

%ORDER: q delta nu

%MCMC particulars
noParam=3;                  %number of parameters for inference
numberIterations=500000;    %number of MCMC iterations
ParamCurrent=log([1e-6,  0.0028, 0.00025]);
noAccept=0;                 %starting acceptance rate
adaptingThresh=2*noParam+1; %adaptive MCMC threshold
chain=ParamCurrent;

%Likelihood for starting conditions
[ seropositive,x0 ] = SimulateMSEIRS( exp(ParamCurrent), Lmax );
%interpolate the seroposivity quartely to match the data
prevalence=interp1(1:76,seropositive,1:0.25:76);

%calculate the binomial log likelihood
LL=0;
parfor index=1:length(age)
    LL=LL+log(binopdf(round(positive(index)),sample_size(index),prevalence(age(index)*4+1)));
end
Like=LL;

%PRIORS
%probability of transmission
q_PRIOR=log(unifpdf(exp(ParamCurrent(1)),0,30));    
%rate of loss of maternal antibodies
delta_PRIOR=log(normpdf(exp(ParamCurrent(2)),0.00275524,0.000430117));    
%rate of loss of immunity
nu_PRIOR=log(normpdf(exp(ParamCurrent(3)),1/(5.1*365),5e-05));

%Acceptance probability
accCurrent=LL+q_PRIOR+delta_PRIOR+nu_PRIOR+ sum(ParamCurrent);  %parameters arelog transformed so add them


%initial proposal distribution
proposal= (0.001^2)*eye(noParam)/noParam;
%begin MCMC
for i=1:numberIterations
    
    %propose values
    if i<adaptingThresh
        
        ParamProp=(mvnrnd(ParamCurrent,proposal));
        
        %%%%MCMCSTEP %%%%
        %LOG LIKELIHOOD
        [ seropositive ] = SimulateMSEIRS( exp(ParamProp), Lmax );
        %interpolate this quarterly
        prevalence=interp1(1:76,seropositive,1:0.25:76);
        
        LL=0;
        for index=1:length(age)
            LL=LL+log(binopdf(round(positive(index)),sample_size(index),prevalence(age(index)*4+1)));
        end
        
        %PRIORS
        q_PRIOR=log(unifpdf(exp(ParamProp(1)),0,30));
        delta_PRIOR=log(normpdf(exp(ParamProp(2)),0.00275524,0.000430117));   
        nu_PRIOR=log(normpdf(exp(ParamProp(3)),1/(5.1*365),5e-05));
        
        %ACCEPTANCE
        accProp=LL+q_PRIOR+delta_PRIOR+nu_PRIOR + sum(ParamProp);
        
        acc=min(1,(exp(accProp-accCurrent)));
        
        
        if rand(1)<acc
            %ACCEPT
            ParamCurrent=ParamProp;
            accCurrent=accProp;
            Like(i+1)=LL;
            noAccept=noAccept+1;
        else
            Like(i+1)=Like(i);
        end
        %%%%%%%%%%%%%%%%%%
        
    else
        
        if rand(1)<0.9 || i==adaptingThresh
            %update covariance matrix
            disp('UPDATE ADAPT')
            adapt=cov(chain(10:end,:))* (2.38 ^ 2) / noParam;
            try
                ParamProp=(mvnrnd(ParamCurrent,adapt));
            catch
                ParamProp=(mvnrnd(ParamCurrent,proposal));
            end
        else
            ParamProp=(mvnrnd(ParamCurrent,proposal));
        end
        
        %%%%MCMCSTEP %%%%
        %LOG LIKELIHOOD
        [ seropositive ] = SimulateMSEIRS( exp(ParamProp), Lmax );
        %interpolate this quarterly
        prevalence=interp1(1:76,seropositive,1:0.25:76);
        %BINOMIAL LOG LIKE
        LL=0;
        for index=1:length(age)
            LL=LL+log(binopdf(round(positive(index)),sample_size(index),prevalence(age(index)*4+1)));
        end
        %PRIORS
        q_PRIOR=log(unifpdf(exp(ParamProp(1)),0,30));
        delta_PRIOR=log(normpdf(exp(ParamProp(2)),0.00275524,0.000430117));  
        nu_PRIOR=log(normpdf(exp(ParamProp(3)),1/(5.1*365),5e-05));
        
        %ACCEPTANCE
        accProp=LL+q_PRIOR+delta_PRIOR+nu_PRIOR + sum(ParamProp);
        
        acc=min(1,(exp(accProp-accCurrent)));
        
        
        if rand(1)<acc
            %ACCEPT
            ParamCurrent=ParamProp;
            accCurrent=accProp;
            Like(i+1)=LL;
            noAccept=noAccept+1;
        else
            Like(i+1)=Like(i);
        end
        %%%%%%%%%%%%%%%%%%
        
    end
    
    chain(i+1,:)=ParamCurrent;
    
    subplot(2,3,1)
    plot(exp(chain(:,1)))
    title('q')
    subplot(2,3,2)
    plot(exp(chain(:,2)))
    title('delta')
    subplot(2,3,3)
    plot(exp(chain(:,3)))
    title('nu')
    subplot(2,3,[4 5 6])
    plot(Like)
    
    drawnow
    disp('acceptance rate is: ')
    disp(noAccept/i);
    
    save('serology_inference')
end
%%

%cut off initial burn in period
chain1=exp(chain(1000:end,:));
LOGLIKE=Like(1000:end);

%plot key features of the chain
figure
scatter(age+1,positive./sample_size,10*sample_size)
hold on
[ seropositive_mode] = SimulateMSEIRS( median(chain1), Lmax);
plot(seropositive_mode,'k');
[ seropositive_low] = SimulateMSEIRS( [quantile(chain1(:,1),0.025),median(chain1(:,[2,3]))], Lmax);
plot(seropositive_low,'k.');
[ seropositive_high] = SimulateMSEIRS( [quantile(chain1(:,1),0.975),median(chain1(:,[2,3]))], Lmax);
plot(seropositive_high,'k.');


%CALCULATING R0
L=120*mean(exp(-cumsum(mu_whole)));     %LIFE EXPECTANCY
q=chain1(:,1);

M=zeros(76);
CSmu=exp(-cumsum(mu));
for i=1:76
    M(i,i)=CSmu(i);
end

for j=1:length(q)
    for i=1:Lmax
        beta(i,:)=q(j)*C(i,:);
    end
    T=[zeros(Lmax),M*beta(:,:)*N/(L);...
        zeros(Lmax),zeros(Lmax)];
    
    
    SIGMA=[alpha*eye(Lmax), zeros(Lmax);...
        -alpha*eye(Lmax) , gamma*eye(Lmax)];
    
    eivalues=eig(T/SIGMA); eivalues=eivalues(eivalues>0);  R0(j)=max(eivalues);
end
figure; histfit(R0,100,'kernel')

%%%trace plot
figure
for i=1:3
    subplot(1,3,i)
    plot(chain1(:,i),'k')
end



%%% priorposteriorplot
figure
%%%q
pd=makedist('uniform','upper',30,'lower',0);
x_values = linspace(min(chain1(:,1)),max(chain1(:,1)));
subplot(1,3,1)
hold on
y = pdf(pd,x_values);
plot(x_values,y,'r','LineWidth',2)

pd=fitdist(chain1(:,1),'kernel');
y = pdf(pd,x_values);
plot(x_values,y,'LineWidth',2)

%%%delta
% pd=makedist('normal','mu',0.00326095,'sigma',0.000388532);  %PELOSI
pd=makedist('normal','mu',0.00275524,'sigma',0.000430117);%SON
x_values = linspace(0,max(chain1(:,2)));
subplot(1,3,2)
hold on
y = pdf(pd,x_values);
plot(x_values,y,'r','LineWidth',2)

pd=fitdist(chain1(:,2),'kernel');
y = pdf(pd,x_values);
plot(x_values,y,'LineWidth',2)

%%%nu
pd=makedist('normal','mu',1/(5.1*365),'sigma',5e-05); %PELOSI
% pd=makedist('normal','mu',0.000255163,'sigma',0.000162234 ); %SON
x_values = linspace(0,1e-3);
subplot(1,3,3)
hold on
y = pdf(pd,x_values);
plot(x_values,y,'r','LineWidth',2)

pd=fitdist(chain1(:,3),'kernel');
y = pdf(pd,x_values);
plot(x_values,y,'LineWidth',2)

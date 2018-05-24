function dydt=T_ode_dynamic(t,Tandh,R,F,inpvec)
    
    
    
    h=Tandh(1)
    T=Tandh(2)
   
    
   %tolerance
options.StepTolerance = 1e-10;
    rho=1000;
    D_T=0.145;
    Ac=0.25*pi()*(D_T)^2;
    delta_H=2100;
    cp=4.1855;
     


    Tc=23;

    hss=0.22;
    Tss=25.7;
    msteamss=18.36/3600;
    mcss=7.12/60;
    
    %input vector
    mc=inpvec(1,1)
    msteam=inpvec(2,1)
    
    
    dydt(1)=-R/(rho*Ac)*h^0.5+1/(rho*Ac)*mc
  
    dydt(2)=(mc*Tc/(rho*Ac)+msteam*delta_H*F/(rho*Ac*cp))*h^-1 -mc/(rho*Ac)*(T*h^-1)

    
    dydt=[dydt(1);dydt(2)]
end
    
function dydt=T_ode_dynamic(t,Tandh,i,j,R1,R2,R3,F)
    
    
    
    h1=Tandh(1)
    h2=Tandh(2)
    h3=Tandh(3)
    T1=Tandh(4)
    T2=Tandh(5)
    T3=Tandh(6)
    
    rho=1000
    D_T=0.145
    Ac=0.25*pi()*(D_T)^2
    cp=4.1855
    %disturbance variables( use S.S. values)
    delta_H=2100
    Tc=4
    Th=50
    
    %input vector
    inpvec=inputs(i,j)
    mh=inpvec(1)
    mc=inpvec(2)
    msteam=inpvec(3)
    
    
    dydt(1)=-R1/(rho*Ac)*h1^0.5+1/(rho*Ac)*mh
    dydt(2)=R1/(rho*Ac)*h1^0.5-R2/(rho*Ac)*h2^0.5
    dydt(3)=-R3/(rho*Ac)*h3^0.5+1/(rho*Ac)*mc
    dydt(4)=(mh*Th/(rho*Ac)+msteam*delta_H*F/(rho*Ac*cp))*h1^-1 -mh/(rho*Ac)*(T1*h1^-1)
    dydt(5)=R1/(rho*Ac)*(T1*h1^0.5*h2^-1) - R1/(rho*Ac)*(T2*h1^0.5*h2^-1)
    dydt(6)=Tc*mc/(rho*Ac)*h3^-1- mc/(rho*Ac)*(T3*h3^-1)
    
    dydt=[dydt(1);dydt(2);dydt(3);dydt(4);dydt(5);dydt(6)]
end
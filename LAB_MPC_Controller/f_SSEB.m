function RHS=f_SSEB(x)
    F=x(1)
    T2ss=x(2)
    T3ss=x(3)
    
    h1ss=0.2
    h2ss=0.2
    h3ss=0.2
    rho=1000
    D_T=0.145
    Ac=0.25*pi()*(D_T)^2
    delta_H=2100
    cp=4.1855
    T1ss=60
    Tc=4
    Th=50
    mh=0.15
    mc=0.1
    msteam=8/3600
    
    R=fsolve(@f_SSMB,[1;1;1])
    R1=R(1)
    R2=R(2)
    R3=R(3)
    
    RHS(1)=(mh*Th+msteam*delta_H*F/(cp))*h1ss^-1 -mh*(T1ss*h1ss^-1)
    RHS(2)=T1ss - T2ss
    RHS(3)=Tc- T3ss
end    
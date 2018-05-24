function f1=f_SSMB(R)
    R1=R(1)
    R2=R(2)
    R3=R(3)
    h1ss=0.2
    h2ss=0.2
    h3ss=0.2
    rho=1000
    D_T=0.145
    Ac=0.25*pi()*(D_T)^2
    mh=0.15
    mc=0.1
    
    f1(1)=-R1*h1ss.^0.5+mh
    f1(2)=R1-R2
    f1(3)=-R3*h3ss.^0.5+mc
    


end
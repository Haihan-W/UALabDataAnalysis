clear all
%% Q2. S.S. h and T
    rho=1000;
    D_T=0.145;
    Ac=0.25*pi()*(D_T)^2;
    delta_H=2100;
    cp=4.1855;
     
    %SS inputs
    mhss=0.15;
    mcss=0.1;
    Thss=50;
    Tcss=4;
    msteamss=8/3600;
    
    %SS states
    h1ss=0.2;
    h2ss=0.2;
    h3ss=0.2;
    T1ss=60;
    
    Output=fsolve(@f_SSEB,[0.5;1;1]);
    %F=Output(1) 
    %%%%%%%% note: using fsolve can't get precies value, so use
    %%%%%%%% explicit form below to obtian F
    F = (mhss*T1ss-mhss*Thss)/(msteamss*delta_H/cp);
    T2ss=Output(2);
    T3ss=Output(3);
    
    %parameters
     %R=fsolve(@f_SSMB,[1;1;1])
    %R1=R(1)
    %R2=R(2)
    %R3=R(3)
    %%%%%%%% note: using fsolve can't get precies value for ^0.5, so use
    %%%%%%%% explicit form below to obtian R1...R3
    R1=mhss/(h1ss.^0.5);
    R2=R1;
    R3=mcss/(h3ss.^0.5);
   
   
 %% Q3. dynamic h and T
    % 2*3 cell (i=1--> +20%; i=2--> -20%)
   Cell = {'mhot+20%', 'mc+20%','msteam+20%'; 'mhot-20%', 'mc-20%','msteam-20%' };
    for i=1:2
        for j=1:3
            
            options = odeset('AbsTol',1e-10,'RelTol',1e-10);
            [t,Tandh]=ode45(@T_ode_dynamic, [0,2000], [h1ss;h2ss;h3ss;T1ss;T2ss;T3ss],options,i,j,R1,R2,R3,F);
                        h1=Tandh(:,1);
                        h2=Tandh(:,2);
                        h3=Tandh(:,3);
                        T1=Tandh(:,4);
                        T2=Tandh(:,5);
                        T3=Tandh(:,6);
                        if i==1
                            figure(j);
                        else
                            figure(3+j);
                        end
                            for count=1:6
                                subplot(2,3,count)
                                if count<=3
                                plot(t,eval(['h' num2str(count) '']))
                                xlabel('t')
                                ylabel(['h' num2str(count) ''])
                                axis([0 2000 0.1 0.3])
                                else 
                                plot(t,eval(['T' num2str(count-3) ''])) 
                                xlabel('t')
                                ylabel(['T' num2str(count-3) ''])
                                    if count==6 %T3
                                    axis([0 2000 0 10])
                                    else
                                    axis([0 2000 55 65])
                                    end
                                end 
                            end
                            suptitle(['h and T vs t @' Cell(i,j) ''])
        end
    end

    %% Q4. linearization
    A=zeros(size(6));
    B=zeros(size(6,3));
        % f1=dh1'/dt=a11*h1'+b11*mh'
        A(1,1)=(-R1/(rho*Ac))*0.5*h1ss^-0.5;
        B(1,1)=(1/(rho*Ac));
        
        %f2=dh2'/dt=a21*h1'+a22*h2'
        A(2,1)=(R1/(rho*Ac))*0.5*h1ss^-0.5;
        A(2,2)=(-R2/(rho*Ac))*0.5*h2ss^-0.5;
        
        %f3=dh3'/dt=a33*h3'+b32*mc'
        A(3,3)= (-R3/(rho*Ac))*0.5*h3ss^-0.5;
        B(3,2)=(1/(rho*Ac));
        
        %f4=dT1'/dt=a41*h1'+a44*T1'+b41*mh'+b43*msteam'
        A(4,1)=(1/(rho*Ac)* mhss*Thss)*(-1)*h1ss^-2+(delta_H*F/(rho*Ac*cp)*msteamss)*(-1)*h1ss^-2+(-1/(rho*Ac)*mhss*T1ss)*(-1)*h1ss^-2;
        A(4,4)=(-1/(rho*Ac)* mhss*h1ss^-1);
        B(4,1)=(1/(rho*Ac)*Thss*h1ss^-1)+(-1/(rho*Ac)*T1ss*h1ss^-1);
        B(4,3)=(delta_H*F/(rho*Ac*cp)*h1ss^-1);
        
        %f5=dT2'/dt=a51*h1'+a52*h2'+a54*T1'+a55*T2'
        A(5,1)=(R1/(rho*Ac)*T1ss*h2ss^-1)*0.5*h1ss^-0.5+(-R1/(rho*Ac)*T2ss*h2ss^-1)*0.5*h1ss^-0.5;
        A(5,2)=(R1/(rho*Ac)*T1ss*h1ss^0.5)*(-1)*h2ss^-2+(-R1/(rho*Ac)*T2ss*h1ss^0.5)*(-1)*h2ss^-2;
        A(5,4)=(R1/(rho*Ac)*h1ss^0.5*h2ss^-1);
        A(5,5)=(-R1/(rho*Ac)*h1ss^0.5*h2ss^-1);
        
        %f6=dT3'/dt=a63*h3'+a66*T3'+b62*mc'
        A(6,3)=(1/(rho*Ac)*Tcss*mcss)*(-1)*h3ss^-2+(-1/(rho*Ac)*T3ss*mcss)*(-1)*h3ss^-2;
        A(6,6)=(-1/(rho*Ac)*mcss*h3ss^-1);
        B(6,2)=(1/(rho*Ac)*Tcss*h3ss^-1)+(-1/(rho*Ac)*T3ss*h3ss^-1);
        
     %% Q5. Discretization
        % sampling time = 1 s
        C=eye(6)
        D=0
        sys=ss(A,B,C,D) 
        sysd = c2d(sys,1) 
        
        Ad=sysd.a;
        Bd=sysd.b;
        Cd=sysd.c;
        Dd=sysd.d;
        
    
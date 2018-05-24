function in=Inputs(i,j)
    mh=0.15
    mc=0.1
    msteam=8/3600
    
    inpvec=[mh,mc,msteam]
    
    if i==1
        inpvec(j)=inpvec(j)*1.2
    end
    if i==2
        inpvec(j)=inpvec(j)*0.8
    end
    in= inpvec
end
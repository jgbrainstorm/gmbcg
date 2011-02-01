pro readconfig,filename,inputdir,catdir,radius,truth

    openr,1,filename
    h=strarr(7)
    readf,1,h
    close,1
inputdir=strmid(h[3],9)
catdir=strmid(h[4],7)
radius=float(strmid(h[5],7))
truth=float(strmid(h[6],6))
return
end

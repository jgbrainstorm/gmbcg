function robust_std,x
    outlier_rem,x,ok,bd,/harsh
    res=stdev(x[ok])
    return,res
end

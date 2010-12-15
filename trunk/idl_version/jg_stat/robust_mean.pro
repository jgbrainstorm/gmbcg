function robust_std,x
    outlier_rem,x,ok,bd,/harsh
    res=mean(x[ok])
    return,res
end

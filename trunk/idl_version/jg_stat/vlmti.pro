function vlmti,z
     ;g=mrdfits('vlmtd_400_imag_v2.10.fit',1)
     ;imag=cspline(g.z,g.imag400,z)
     ;the new one is based on mock v3.04
     coef = [ 7.24235186, -29.70849265,  47.82343782, -33.93019772,4.29971516,  10.31565436,  16.54869522]
     zmag=poly(z,coef)
     return,zmag
end

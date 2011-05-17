function vlmti,z
     g=mrdfits('vlmtd_400_imag_v2.10.fit',1)
     imag=cspline(g.z,g.imag400,z)
     return,imag
end

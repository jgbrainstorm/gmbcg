pro ridgeline_evolution

!p.multi=[0,2,2]

x=where(b.gr_ridge eq 1)
y=where(b.ri_ridge eq 1)
z=where(b.iz_ridge eq 1)

x=where(b.photoz le 0.35)
y=where(b.photoz ge 0.35 and b.photoz le 0.7)
z=where(b.photoz ge 0.7 and b.photoz le 1.0)



plotbin,b[x].photoz,b[x].gm_gmr,bin=0.03,/linefit,xtitle='photoz',ytitle='g-r ridgeline',charsize=2
plotbin,b[x].photoz,b[x].gmr,bin=0.03,color=!red,/ov
plotbin,b[y].photoz,b[y].gm_rmi,bin=0.03,/linefit,xtitle='photoz',ytitle='r-i ridgeline',charsize=2
plotbin,b[y].photoz,b[y].rmi,bin=0.03,color=!red,/ov
plotbin,b[z].photoz,b[z].gm_imz,bin=0.03,/linefit,xtitle='photoz',ytitle='i-z ridgeline',charsize=2
plotbin,b[z].photoz,b[z].imz,bin=0.03,color=!red,/ov






plotbin,b[x].photoz,b[x].gm_gmr_wdh,bin=0.05,/linefit,xtitle='photoz',ytitle='g-r ridgeline width',charsize=2
plotbin,b[y].photoz,b[y].gm_rmi_wdh,bin=0.05,/linefit,xtitle='photoz',ytitle='r-i ridgeline width',charsize=2
plotbin,b[z].photoz,b[z].gm_imz_wdh,bin=0.05,/linefit,xtitle='photoz',ytitle='i-z ridgeline width',charsize=2






end

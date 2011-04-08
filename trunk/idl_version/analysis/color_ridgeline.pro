function gmmridgeplot,bcg,z1,z2
         
         x=where(bcg.photoz ge z1 and bcg.photoz le z2)

         if z1 ge 0.1 and z2 le 0.35 then begin
             ok=where(bcg[x].gm_gmr ne 0)
             x=x[ok]
             alpha=[median(bcg[x].gm_mix_gmr_clr),median(bcg[x].gm_mix_gmr_bgd)]
             mu=[mean(bcg[x].gm_gmr),mean(bcg[x].gm_gmr_bgd)]
             sigma=[mean(bcg[x].gm_gmr_wdh),mean(bcg[x].gm_gmr_wdh_bgd)]
             gmmplot,alpha,mu,sigma,xrange=[-0.5,2],npoint=20000,xtitle='g-r',charsize=2
             legend,['z: '+ntostr(z1,5)+' - '+ntostr(z2,5)],box=0,charsize=2
         endif 


         if z1 ge 0.35 and z2 le 0.7 then begin
             ok=where(bcg[x].gm_rmi ne 0)
             x=x[ok]
             alpha=[median(bcg[x].gm_mix_rmi_clr),median(bcg[x].gm_mix_rmi_bgd)]
             mu=[mean(bcg[x].gm_rmi),median(bcg[x].gm_rmi_bgd)]
             sigma=[mean(bcg[x].gm_rmi_wdh),mean(bcg[x].gm_rmi_wdh_bgd)]
             gmmplot,alpha,mu,sigma,xrange=[-0.5,2],npoint=20000,xtitle='r-i',charsize=2
             legend,['z: '+ntostr(z1,5)+' - '+ntostr(z2,5)],box=0,charsize=2
         endif 

         if z1 ge 0.7 and z2 le 1.0 then begin
             ok=where(bcg[x].gm_imz ne 0)
             x=x[ok]
             alpha=[median(bcg[x].gm_mix_imz_clr),median(bcg[x].gm_mix_imz_bgd)]
             mu=[mean(bcg[x].gm_imz),mean(bcg[x].gm_imz_bgd)]
             sigma=[mean(bcg[x].gm_imz_wdh),mean(bcg[x].gm_imz_wdh_bgd)]
             gmmplot,alpha,mu,sigma,xrange=[-0.5,2],npoint=20000,xtitle='i-z',charsize=2
             legend,['z: '+ntostr(z1,5)+' - '+ntostr(z2,5)],box=0,charsize=2
         endif 

         return,0
end



pro color_ridgeline,bcg

   ; bcg=mrdfits('/data/des_mock_catalog/v2.00/gmbcg_cluster/true_bcg/des_mock_v2.00_BCG_gmbcg_v2.5.fit',1)
   ; bcg=mrdfits('/data/des_mock_catalog/v2.10/gmbcg_cluster/vmaglim_i/des_mock_v2.10_BCG_blended_gmbcg_v2.5.fit',1)
   ; bcg=bcg[where(bcg.gm_ngals_weighted ge 10)]
   ; bcg=mrdfits('/data/des_mock_catalog/v2.10/gmbcg_cluster/truebcg_z/vmaglimi/des_mock_v2.10_BCG_blended_gmbcg_v2.5.fit',1)
   ; bcg=mrdfits('/data/des_dc5/DC5_ra_le_335_dec_le_minus_45_ridgeline.fit',1)
   ; bcg=mrdfits('/data/des_dc5/DC5B_ridgeline.fit',1)
   ; bcg=mrdfits('/data/des_mock_catalog/v2.11/gmbcg_cluster/arborz/des_mock_v2.10_BCG_blended_gmbcg_v2.5.fit',1)

  ; bcg=mrdfits('/data/des_dc5/DC5B2/gmbcg_cluster/galaxy/des_mock_v2.10_BCG_blended_gmbcg_v2.5.fit',1)
   bcg=mrdfits('/data/des_dc5/DC5B2/gmbcg_cluster/truth/des_mock_v2.10_BCG_blended_gmbcg_v2.5.fit',1)
 
    bcg=bcg[where(bcg.gm_nn eq 2 and bcg.ngals ge 8)]
    !p.multi=[0,3,4]
    window,xsize=1000,ysize=1100

    t=gmmridgeplot(bcg,0.1,0.15)
    t=gmmridgeplot(bcg,0.15,0.17)
    t=gmmridgeplot(bcg,0.17,0.21)
    t=gmmridgeplot(bcg,0.21,0.24)
    t=gmmridgeplot(bcg,0.24,0.27)
    t=gmmridgeplot(bcg,0.27,0.30)
    t=gmmridgeplot(bcg,0.30,0.33)
    t=gmmridgeplot(bcg,0.33,0.35)
    t=gmmridgeplot(bcg,0.35,0.38)
    t=gmmridgeplot(bcg,0.38,0.42)
    t=gmmridgeplot(bcg,0.42,0.46)
    t=gmmridgeplot(bcg,0.46,0.50)

     
    t=gmmridgeplot(bcg,0.50,0.54)
    t=gmmridgeplot(bcg,0.54,0.60)
    t=gmmridgeplot(bcg,0.60,0.64)
    t=gmmridgeplot(bcg,0.64,0.70)
    t=gmmridgeplot(bcg,0.70,0.74)
    t=gmmridgeplot(bcg,0.74,0.80)
    t=gmmridgeplot(bcg,0.80,0.84)
    t=gmmridgeplot(bcg,0.84,0.90)
    t=gmmridgeplot(bcg,0.90,0.94)
    t=gmmridgeplot(bcg,0.94,0.97)
    t=gmmridgeplot(bcg,0.97,1.)




end

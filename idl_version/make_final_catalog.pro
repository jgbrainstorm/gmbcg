;this pro choose the final catalog for comparison project. 
; the format is per
; https://sites.google.com/site/descwg/Home/cluster-comparison/2009dec-meeting

pro make_final_catalog,catdir,version

   ; cat_dir='/data/des_mock_catalog/v2.10/gmbcg_cluster/vmaglim_i/'  

    bg=mrdfits(catdir+'des_mock_v'+ntostr(version,4)+'_BCG_blended_gmbcg_v2.5.fit',1)
     


    t1=create_struct('rank',0L,'ra',0.D,'dec',0.D,'z',0.,'photoz',0.,'ngals',0,'nfw_lh',0.,'gm_ngals_weighted',0.,'objid',0L,'bic1',0.,'bic2',0.)
    t1=replicate(t1,n_elements(bg.(0)))
    rank=lindgen(n_elements(bg.(0)))
    s=reverse(sort(bg.ngals))
    bg=bg[s]
    t1.rank=rank  
    t1.ra=bg.ra
    t1.dec=bg.dec
    t1.z=bg.photoz
    t1.photoz=bg.photoz
    t1.ngals=bg.ngals
    t1.nfw_lh=bg.nfw_lh
    t1.gm_ngals_weighted=bg.gm_ngals_weighted
    t1.objid = bg.objid
    t1.bic1 = bg.bic1
    t1.bic2 = bg.bic2
    mwrfits,t1,catdir+'des_mock_v'+ntostr(version,4)+'_gmbcg_Hao.fit',/create



end


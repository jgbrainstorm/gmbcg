;This loop call the cluster finder main engine. 

pro gmbcg_loop,input_dir,cat_dir,radius,truth=truth

    ;input_dir='/data/des_mock_catalog/v2.13/original_cat/'
    ;cat_dir='/data/des_mock_catalog/v2.13/gmbcg_cluster/newsvn/'
    
   ; readconfig,filename,input_dir,cat_dir,srad

   ; input_dir='/archive/staging/DES/SIMULATIONS/wechsler/v2.13/'
   ; cat_dir='/archive/staging/DES/SIMULATIONS/wechsler/v2.13/'
    version = 1.0
    if keyword_set(truth) then file=findfile(input_dir+'DES_Mock_v2.13_Baseline_truth*.fit') else file=findfile(input_dir+'DES_Mock_v2.13_Baseline_0*.fit')  
    for patch=0,1 do begin
        print,patch
        
        gal=mrdfits(file[patch],1)
        num=n_elements(gal.(0))
        tt=create_struct('objid',long64(0),'ra',0.D,'dec',0.D,'amag',[0.,0.,0.,0.,0.],'tmag',[0.,0.,0.,0.,0.],'omag',[0.,0.,0.,0.,0.],'omag_err',[0.,0.,0.,0.,0.],'z',0.,'photoz',0.,'photoz_err',0.,'gmr',0.,'rmi',0.,'imz',0.,'zmy',0.,'gmr_err',0.,'rmi_err',0.,'imz_err',0.,'zmy_err',0.,'gm_gmr',0.,'gm_gmr_wdh',0.,'gm_rmi',0.,'gm_rmi_wdh',0.,'gm_imz',0.,'gm_imz_wdh',0.,'gm_zmy',0.,'gm_zmy_wdh',0.,'gr_ridge',0,'ri_ridge',0,'iz_ridge',0,'zy_ridge',0,'nfw_lh',0.,'bcgmag_lh',0.,'used',0,'rcenter',0.,'bcg_gr_lh',0.,'bcg_ri_lh',0.,'bcg_iz_lh',0.,'lh',0.,'bcglh',0.,'GM_NN',0,'ngals',0,'ngals_r200',0,'central',0,'arborz',0.,'arborz_err',0.,'annz',0.,'annz_err',0.,'photoz_gaussian',0.,'gm_mix_gmr_clr',0.,'gm_mix_gmr_bgd',0.,'gm_mix_rmi_clr',0.,'gm_mix_rmi_bgd',0.,'gm_mix_imz_clr',0.,'gm_mix_imz_bgd',0.,'gm_mix_zmy_clr',0.,'gm_mix_zmy_bgd',0.,'isbcg',0.,'GM_gmr_bgd',0.,'GM_gmr_wdh_bgd',0.,'GM_rmi_bgd',0.,'GM_rmi_wdh_bgd',0.,'GM_imz_bgd',0.,'GM_imz_wdh_bgd',0.,'GM_zmy_bgd',0.,'GM_zmy_wdh_bgd',0.,'lim_i',0.,'Ntot',0.,'GM_Ngals_weighted',0.)
        str=replicate(tt,n_elements(gal.(0)))
        str.objid = gal.id
        if keyword_set(truth) then begin
           str.amag = gal.amag
           str.tmag = gal.tmag
           str.z=gal.z
           str.photoz=gal.z
           str.photoz_err=0.
           str.central=gal.central
        endif else begin
            str.omag[0] = gal.mag_g
            str.omag[1] = gal.mag_r
            str.omag[2] = gal.mag_i
            str.omag[3] = gal.mag_z
            str.omag[4] = gal.mag_y
       
            str.omag_err[0] = gal.magerr_g
            str.omag_err[1] = gal.magerr_r
            str.omag_err[2] = gal.magerr_i
            str.omag_err[3] = gal.magerr_z
            str.omag_err[4] = gal.magerr_y
            str.photoz =gal.annz
            str.photoz_err = gal.annz_err
        endelse

        str.ra = double(gal.ra)
        str.dec = double(gal.dec)        
        str.gmr=str.omag[0]-str.omag[1]
        str.rmi=str.omag[1]-str.omag[2]
        str.imz=str.omag[2]-str.omag[3]
        str.zmy=str.omag[3]-str.omag[4]
        str.gmr_err=sqrt(str.omag_err[0]^2+str.omag_err[1]^2)
        str.rmi_err=sqrt(str.omag_err[1]^2+str.omag_err[2]^2)
        str.imz_err=sqrt(str.omag_err[2]^2+str.omag_err[3]^2)
        str.zmy_err=sqrt(str.omag_err[3]^2+str.omag_err[4]^2)
        x_gr=where(str.photoz le 0.35)
        x_ri=where(str.photoz gt 0.35 and str.photoz le 0.7)
        x_iz=where(str.photoz gt 0.7 and str.photoz le 1.0)
        x_zy=where(str.photoz gt 1.0)
        str[x_gr].gr_ridge=1
        str[x_ri].ri_ridge=1
        str[x_iz].iz_ridge=1
        str[x_zy].zy_ridge=1
       ; str.lim_i=vlmti(str.photoz)
        str.lim_i=limi(str.photoz)
        gal=str[where(str.photoz ge 0.1 and str.photoz le 1.0)]
            
        des_mock_gmbcg,cat_dir,gal,radius,patch,version
      
    
    endfor

    gmbcg_des_combine_cat,cat_dir,version  ; combine the BCG catalogs
    des_gmbcg_percolation,cat_dir,version   
    make_final_catalog, cat_dir,version
    gmbcg_des_member_combine,cat_dir,version
  

end

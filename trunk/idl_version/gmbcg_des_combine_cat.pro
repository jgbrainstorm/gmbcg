;only combine the bcg catalog

pro gmbcg_des_combine_cat,cat_dir,version

     bgf=findfile(cat_dir+'des_mock_v'+ntostr(version,4)+'_BCG_gmbcg*')
     bg=mrdfits(bgf[0],1)
     for i=1,n_elements(bgf)-1 do begin
         b=mrdfits(bgf[i],1)
         bg=[bg,b] 
     endfor
   
     mwrfits,bg,cat_dir+'des_mock_v'+ntostr(version,4)+'_BCG_gmbcg_v2.5.fit',/create

end

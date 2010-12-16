pro gmbcg_des_member_combine,cat_dir,version

    bg=mrdfits(cat_dir+'des_mock_v'+ntostr(version,4)+'_gmbcg_Hao.fit',1)
    
    bmf=findfile(cat_dir+'des_mock_v'+ntostr(version,4)+'_BCGMB_gmbcg*')
    
    tp=create_struct('rank',0L,'id',0L)
    tp=replicate(tp,n_elements(bg.(0)))
    tp.rank=bg.rank
    tp.id=long(bg.objid)
    
    for i=1,n_elements(bmf)-1 do begin

        bm=mrdfits(bmf[i],1)
        add_tags,bm,['rank'],['0L'],bmb
        bm=0
        close_match_radec,bg.ra,bg.dec,bmb.bcgra,bmb.bcgdec,i1,i2,0.0004,1
        bmb[i2].rank=bg[i1].rank
        nm=n_elements(i2)
        tmp=create_struct('rank',0L,'id',0L)
        tmp=replicate(tmp,nm)
        tmp.rank=bmb[i2].rank
        tmp.id=long(bmb[i2].objid)
        tp=[tp,tmp]

    endfor
   
    mwrfits,tp,cat_dir+'des_mock_v'+ntostr(version,4)+'_gmbcg_Hao.fit'

end
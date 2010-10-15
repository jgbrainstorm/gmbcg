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
        t=inner_join(bg.objid,bmb.bcgid)
        bmb[t.match2].rank=bg[t.match1].rank
        nm=n_elements(t.match2)
        tmp=create_struct('rank',0L,'id',0L)
        tmp=replicate(tmp,nm)
        tmp.rank=bmb[t.match2].rank
        tmp.id=long(bmb[t.match2].objid)
        tp=[tp,tmp]

    endfor
   
    mwrfits,tp,cat_dir+'des_mock_v'+ntostr(version,4)+'_gmbcg_Hao.fit'

end

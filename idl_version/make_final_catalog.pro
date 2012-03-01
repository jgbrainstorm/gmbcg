;this pro choose the final catalog for comparison project. 
; the format is per
; https://sites.google.com/site/descwg/Home/cluster-comparison/2009dec-meeting

function rerank,bcg,z1,z2,n
    idx = where(bcg.photoz ge z1 and bcg.photoz lt z2)
    s = reverse(sort(bcg[idx].ngals))
    sortidx = idx[s[0:n-1]]
    return, sortidx
end

pro make_final_catalog,catdir,version

    bg=mrdfits(catdir+'des_mock_v'+ntostr(version,4)+'_BCG_gmbcg_v2.5.fit',1)

    t1=create_struct('rank',0L,'ra',0.D,'dec',0.D,'z',0.,'photoz',0.,'ngals',0,'nfw_lh',0.,'gm_ngals_weighted',0.,'objid',0L,'bic1',0.,'bic2',0.)
    t1=replicate(t1,n_elements(bg.(0)))
    t1.rank=-999  
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

    num = 260
    idx = rerank(t1,0.05,0.25,num)
    rank = lindgen(num)+1
    t1[idx].rank = rank

    num = 531
    idx = rerank(t1,0.25,0.4,num)
    rank = lindgen(num)+1
    t1[idx].rank = rank

    num = 1194
    idx = rerank(t1,0.4,0.6,num)
    rank = lindgen(num)+1
    t1[idx].rank = rank

    
    num = 1091
    idx = rerank(t1,0.6,0.75,num)
    rank = lindgen(num)+1
    t1[idx].rank = rank

    num = 1532
    idx = rerank(t1,0.75,1.,num)
    rank = lindgen(num)+1
    t1[idx].rank = rank

    num = 686
    idx = rerank(t1,1.,1.15,num)
    rank = lindgen(num)+1
    t1[idx].rank = rank

    num = 790
    idx = rerank(t1,1.15,1.4,num)
    rank = lindgen(num)+1
    t1[idx].rank = rank

    t1=t1[where(t1.rank ne -999)]
    mwrfits,t1,catdir+'des_mock_v'+ntostr(version,4)+'_gmbcg_Hao.fit',/create

end


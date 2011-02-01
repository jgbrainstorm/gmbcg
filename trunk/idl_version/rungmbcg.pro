pro run_gmbcg

    readconfig,'config',input_dir,cat_dir,radius,truth
    if (truth eq 0) then begin
        gmbcg_loop,input_dir,cat_dir,radius
    endif else begin
        gmbcg_loop,input_dir,cat_dir,radius,/truth
    endelse
    return 
end

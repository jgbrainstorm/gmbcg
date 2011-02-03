pro rungmbcginner

    readconfig,'config',input_dir,cat_dir,radius,tru
    if (tru eq 0) then begin
        gmbcg_loop,input_dir,cat_dir,radius
    endif else begin
        gmbcg_loop,input_dir,cat_dir,radius,/truth
    endelse
    return 
end

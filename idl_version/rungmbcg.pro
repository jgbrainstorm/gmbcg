pro run_gmbcg

    readconfig,'config',input_dir,cat_dir,radius
    gmbcg_loop,input_dir,cat_dir,radius
    ;if you want to run on the truth table, use:
    ;gmbcg_loop,input_dir,cat_dir,radius,/truth
    return 
end

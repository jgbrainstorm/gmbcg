pro select_des_mock_red,b,ind

ind=where(b.imz ge 0.1 and b.imz ge b.rmi*0.5-0.3 and b.imz le b.rmi*2.08-0.08 and b.imz le b.rmi*(-0.3)+1.7 and b.gmr ge 0.1 and b.gmr le 4 and b.rmi le 0.91*b.gmr+0.82 and b.rmi ge 0.1 and b.rmi le 2 and b.rmi ge 0.4*b.gmr-0.3 and b.rmi le -0.4*b.gmr+3 and b.zmy ge 0.214*b.imz-0.15 and b.zmy le 0.7*b.imz+0.1 and b.zmy le (-0.5)*b.imz +1.3)


return

end

set(SOURCES 	d1mach.c 
	da7sst.c dc7vfn.c dd7dgb.c dd7dog.c dd7dup.c dd7mlp.c 
	dd7tpr.c dd7upd.c df7dhb.c df7hes.c dg7itb.c dg7lit.c 
	dg7qsb.c dg7qts.c dh2rfa.c dh2rfg.c ditsum.c divset.c 
	dl7itv.c dl7ivm.c dl7msb.c dl7mst.c dl7nvr.c dl7sqr.c 
	dl7srt.c dl7svn.c dl7svx.c dl7tsq.c dl7tvm.c dl7upd.c 
	dl7vml.c dmngb.c  dmnhb.c  dn2cvp.c dn2gb.c  dn2g.c 
	dn2lrd.c dn2rdp.c do7prd.c dparck.c dq7apl.c dq7rad.c 
	dq7rgs.c dq7rsh.c dr7mdc.c dr7tvm.c drldst.c drmngb.c 
	drmnhb.c drn2gb.c drn2g.c  ds7bqn.c ds7dmp.c ds7ipr.c 
	ds7lup.c ds7lvm.c dv2axy.c dv2nrm.c dv7cpy.c dv7dfl.c 
	dv7ipr.c dv7scl.c dv7scp.c dv7shf.c dv7swp.c dv7vmp.c 
	dw7zbf.c i1mach.c i7copy.c i7mdcn.c i7pnvr.c i7shft.c 
	stopx.c)

	add_library(cport ${SOURCES})
	target_link_libraries(cport PUBLIC f2c)
	if(MSVC)
		target_compile_options(cport PUBLIC /wd4554)
	endif()
jlib_shared_mkl:
	icc -fPIC -shared -o libjulia.so mkl_.cpp -qmkl

docker:
	icc -fPIC -shared -o libmkl.julia.so mkl_.cpp -qmkl

jlib_static_mkl:
	icc -fPIC -shared -o libjulia.so mkl_.cpp -L/opt/intel/oneapi/mkl/2023.1.0/lib/intel64/libmkl_intel_lp64.a -Wl,--start-group /opt/intel/oneapi/mkl/2023.1.0/lib/intel64/libmkl_intel_lp64.a /opt/intel/oneapi/mkl/2023.1.0/lib/intel64/libmkl_sequential.a /opt/intel/oneapi/mkl/2023.1.0/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl

windows:
	cl /EHsc /LD mkl_.cpp /link /LIBPATH:"C:\Program Files (x86)\Intel\oneAPI\mkl\2022.0.0\lib" mkl_rt.lib /OUT:libjulia.dll

cleanj:
	rm libjulia.so

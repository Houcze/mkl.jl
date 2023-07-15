jlib_shared_mkl:
	icc -fPIC -shared -o libjulia.so mkl_.cpp -qmkl

jlib_static_mkl:
	icc -fPIC -shared -o libjulia.so mkl_.cpp -L/opt/intel/oneapi/mkl/2023.1.0/lib/intel64/libmkl_intel_lp64.a -Wl,--start-group /opt/intel/oneapi/mkl/2023.1.0/lib/intel64/libmkl_intel_lp64.a /opt/intel/oneapi/mkl/2023.1.0/lib/intel64/libmkl_sequential.a /opt/intel/oneapi/mkl/2023.1.0/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl

cleanj:
	rm libjulia.so

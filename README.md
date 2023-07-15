# mkl.jl
Intel MKL for Julia (build for homework usage)

* Usage
Download this codebase to the current path. If you're using Ubuntu (as of July 15, 2023, the Ubuntu version on Colab), you should rename `lib.ubuntu2004` to `lib`:
```sh
mv lib.ubuntu2004 lib
```
Run the following command:
```sh
. set ld_library_path.sh
```
This will add `lib` to the `LD_LIBRARY_PATH`. 

Next, run:
```sh
julia test_mkl.jl
```
to check if you can use the precompiled dynamic library. 

If the installation is successful, you will see the following output:
```
lambda 1 is -2.0
corresponding eigenvector: [-0.7593479889713342, -0.34457275170805207, 0.5519603703396933]
lambda 2 is -1.9999999999999982
corresponding eigenvector: [-0.043027173021516596, 0.9020385411581769, -0.42950568406833034]
lambda 3 is 4.000000000000002
corresponding eigenvector: [-0.5773502691896258, 0.5773502691896258, -0.5773502691896258]
[1.0 -3.0 3.0; 3.0 -5.0 3.0; 6.0 -6.0 4.0]
[2.0 6.0 -8.0; 12.0 1.0 5.0; -16.0 -43.0 3.0]
[1.0 0.0 0.0; 0.0 1.0 0.0; -0.0 -0.0 1.0]
```
Now you can use `mkl.jl`.

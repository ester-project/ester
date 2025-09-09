{
  lib,
  stdenv,
  cmake,
  pkg-config,
  swig,
  binutils,
  which,
  gfortran,
  blas,
  lapack,
  hdf5-cpp,
  python3,
  llvmPackages,
  makeBinaryWrapper,
}: let
  pythonWithPackages = python3.withPackages (p: [p.numpy p.matplotlib p.tkinter]);
  pythonPath = "$out/lib/python:${pythonWithPackages}/${pythonWithPackages.sitePackages}";
  pythonExecutable = "${pythonWithPackages}/bin/${pythonWithPackages.executable}";
in
  stdenv.mkDerivation {
    pname = "ester";
    version = "0";

    src = lib.fileset.toSource {
      root = ./.;
      fileset = lib.fileset.gitTracked ./.;
    };

    cmakeFlags = [
      "-DPYTHON_SITE=lib/python" 
      "-DCBLAS_LIBRARIES=cblas"
      "-DCMAKE_INSTALL_PREFIX=${placeholder "out"}"
    ];
    preConfigure = ''
      numpy_include=$(${pythonExecutable} -c "import numpy; print(numpy.get_include())")
      cmakeFlagsArray+=("-DPYTHON_NUMPY_INCLUDE_DIR=$numpy_include")
    '';

    dontStrip = stdenv.isDarwin;

    nativeBuildInputs = [
      cmake
      pkg-config
      swig
      binutils
      which
      gfortran
      makeBinaryWrapper
    ];

    buildInputs =
      [
        blas
        lapack
        hdf5-cpp.dev
        pythonWithPackages
      ]
      ++ lib.optional stdenv.cc.isClang llvmPackages.openmp;

    #postPatch = ''
    #  substituteInPlace python/CMakeLists.txt --replace 'execute_process(' 'set(PYTHON_SITE "$out/lib/python")\n# execute_process('
    #'';

    postInstall =
      /*
      bash
      */
      ''
        for prog in ester_info gen_output star1d star2d star_evol version; do
          wrapProgram  "$out/bin/$prog" \
            --set NIX_PYTHONPREFIX "${python3}" \
            --set NIX_PYTHONEXECUTABLE ${pythonExecutable} \
            --set NIX_PYTHONPATH ${pythonPath} \
            --set PYTHONNOUSERSITE "true"
        done
      '';
  }

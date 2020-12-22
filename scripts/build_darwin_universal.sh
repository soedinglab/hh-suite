#!/bin/sh -ex
abspath() {
    if [ -d "$1" ]; then
        (cd "$1"; pwd)
    elif [ -f "$1" ]; then
        if [ -z "${1##*/*}" ]; then
            echo "$(cd "${1%/*}"; pwd)/${1##*/}"
        else
            echo "$(pwd)/$1"
        fi
    elif [ -d "$(dirname "$1")" ]; then
        echo "$(cd "$(dirname "$1")"; pwd)/$(basename "$1")"
    fi
}

REPO="$(abspath "$1")"
BUILD="$(abspath "$2")"
CPUS=${3:-$(nproc --all 2>/dev/null || sysctl -n hw.ncpu)}

if [ ! -d "$REPO" ]; then
    echo "${$REPO} repository missing"
    exit 1
fi

export MACOSX_DEPLOYMENT_TARGET=10.13

mkdir -p "$BUILD/build_libomp" && cd "$BUILD/build_libomp"
wget -q http://github.com/llvm/llvm-project/releases/download/llvmorg-11.0.0/openmp-11.0.0.src.tar.xz
tar xvf openmp-11.0.0.src.tar.xz
cd openmp-11.0.0.src
wget https://raw.githubusercontent.com/Homebrew/formula-patches/7e2ee1d7/libomp/arm.patch
patch -p1 < arm.patch

mkdir -p "$BUILD/build_libomp/openmp-11.0.0.src/build-amd64" && cd "$BUILD/build_libomp/openmp-11.0.0.src/build-amd64"
cmake -DLIBOMP_ENABLE_SHARED=OFF -DLIBOMP_INSTALL_ALIASES=OFF -DLIBOMP_ARCH=x86_64 -DCMAKE_CXX_FLAGS="-arch x86_64" ..
make -j${CPUS}
export LIBOMP_AMD64="$BUILD/build_libomp/openmp-11.0.0.src/build-amd64/runtime/src"

mkdir -p "$BUILD/build_sse41" && cd "$BUILD/build_sse41"
cmake -DCMAKE_INSTALL_PREFIX="${BUILD}/hhsuite" -DCMAKE_BUILD_TYPE=Release -DCHECK_MPI=0 -DHAVE_SSE4_1=1 -DCMAKE_C_FLAGS="-arch x86_64" -DCMAKE_CXX_FLAGS="-arch x86_64" -DBUILD_SHARED_LIBS=OFF -DCMAKE_FIND_LIBRARY_SUFFIXES=".a" -DOpenMP_C_FLAGS="-Xpreprocessor -fopenmp -I${LIBOMP_AMD64}" -DOpenMP_C_LIB_NAMES=omp -DOpenMP_CXX_FLAGS="-Xpreprocessor -fopenmp -I${LIBOMP_AMD64}" -DOpenMP_CXX_LIB_NAMES=omp -DOpenMP_omp_LIBRARY=${LIBOMP_AMD64}/libomp.a "$REPO"
make -j${CPUS}
make install

mkdir -p "$BUILD/build_avx2" && cd "$BUILD/build_avx2"
cmake -DCMAKE_BUILD_TYPE=Release -DCHECK_MPI=0 -DHAVE_AVX2=1 -DCMAKE_C_FLAGS="-arch x86_64h" -DCMAKE_CXX_FLAGS="-arch x86_64h" -DBUILD_SHARED_LIBS=OFF -DCMAKE_FIND_LIBRARY_SUFFIXES=".a" -DOpenMP_C_FLAGS="-Xpreprocessor -fopenmp -I${LIBOMP_AMD64}" -DOpenMP_C_LIB_NAMES=omp -DOpenMP_CXX_FLAGS="-Xpreprocessor -fopenmp -I${LIBOMP_AMD64}" -DOpenMP_CXX_LIB_NAMES=omp -DOpenMP_omp_LIBRARY=${LIBOMP_AMD64}/libomp.a "$REPO"
make -j${CPUS}

export MACOSX_DEPLOYMENT_TARGET=11.0

mkdir -p "$BUILD/build_libomp/openmp-11.0.0.src/build-arm64" && cd "$BUILD/build_libomp/openmp-11.0.0.src/build-arm64"
cmake -DLIBOMP_ENABLE_SHARED=OFF -DLIBOMP_INSTALL_ALIASES=OFF -DLIBOMP_ARCH=aarch64 -DCMAKE_CXX_FLAGS="-arch arm64" -DLIBOMP_ASMFLAGS="-arch arm64" ..
make -j${CPUS}
export LIBOMP_AARCH64="$BUILD/build_libomp/openmp-11.0.0.src/build-arm64/runtime/src"

mkdir -p "$BUILD/build_arm64" && cd "$BUILD/build_arm64"
cmake -DCMAKE_BUILD_TYPE=Release \
      -DCHECK_MPI=0 \
      -DHAVE_ARM8=1 \
      -DCMAKE_C_FLAGS="-arch arm64" \
      -DCMAKE_CXX_FLAGS="-arch arm64" \
      -DBUILD_SHARED_LIBS=OFF \
      -DCMAKE_FIND_LIBRARY_SUFFIXES=".a" \
      -DOpenMP_C_FLAGS="-Xpreprocessor -fopenmp -I${LIBOMP_AARCH64}" \
      -DOpenMP_C_LIB_NAMES=omp \
      -DOpenMP_CXX_FLAGS="-Xpreprocessor -fopenmp -I${LIBOMP_AARCH64}" \
      -DOpenMP_CXX_LIB_NAMES=omp \
      -DOpenMP_omp_LIBRARY=${LIBOMP_AARCH64}/libomp.a \
      "$REPO"
make -j${CPUS}

for BINARY in "${BUILD}/hhsuite/bin/"*; do
    BINARY_NAME="$(basename "${BINARY}")"
    case "${BINARY_NAME}" in
        *ffindex*)
          lipo -create \
              -arch x86_64 "$BUILD/build_sse41/lib/ffindex/src/${BINARY_NAME}" \
              -arch x86_64h "$BUILD/build_avx2/lib/ffindex/src/${BINARY_NAME}" \
              -arch arm64 "$BUILD/build_arm64/lib/ffindex/src/${BINARY_NAME}" \
              -output "$BUILD/hhsuite/bin/${BINARY_NAME}"
          ;;
        *)
          lipo -create \
              -arch x86_64 "$BUILD/build_sse41/src/${BINARY_NAME}" \
              -arch x86_64h "$BUILD/build_avx2/src/${BINARY_NAME}" \
              -arch arm64 "$BUILD/build_arm64/src/${BINARY_NAME}" \
              -output "$BUILD/hhsuite/bin/${BINARY_NAME}"
          ;;
    esac
done



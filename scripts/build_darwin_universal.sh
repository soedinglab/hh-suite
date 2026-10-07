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
    echo "${REPO} repository missing"
    exit 1
fi

ALLOWED_DL_LIBS="lib(System\.B|c\+\+|objc)\."

check_binary() {
    otool -L "$1"
    if [ "$(otool -L "$1" | tail -n +2 | grep -v -E "${ALLOWED_DL_LIBS}" )" != "" ]; then
        echo "Too many linked libraries found in $1 binary. Build is not static!"
        exit 1
    fi

    if ! vtool -show "$1" | tee | grep minos | \
         awk -v version="${MACOSX_DEPLOYMENT_TARGET}" '$2 > version { exit 1 }'
    then
        echo "macOS deployment target was not set correctly"
        exit 1
    fi
}

export MACOSX_DEPLOYMENT_TARGET=10.15

mkdir -p "$BUILD/build_libomp" && cd "$BUILD/build_libomp"
OMPVERSION=20.1.7
wget -qO- https://github.com/llvm/llvm-project/releases/download/llvmorg-${OMPVERSION}/cmake-${OMPVERSION}.src.tar.xz | tar xJvf -
wget -qO- https://github.com/llvm/llvm-project/releases/download/llvmorg-${OMPVERSION}/openmp-${OMPVERSION}.src.tar.xz | tar xJvf -
mv cmake-${OMPVERSION}.src cmake

mkdir -p "$BUILD/build_libomp/openmp-${OMPVERSION}.src/build-amd64" && cd "$BUILD/build_libomp/openmp-${OMPVERSION}.src/build-amd64"
cmake \
    -DLIBOMP_ENABLE_SHARED=OFF \
    -DLIBOMP_INSTALL_ALIASES=OFF \
    -DLIBOMP_ARCH=x86_64 \
    -DCMAKE_OSX_ARCHITECTURES=x86_64 \
    -DCMAKE_C_FLAGS="-arch x86_64" \
    -DCMAKE_CXX_FLAGS="-arch x86_64" \
    -DLIBOMP_ASMFLAGS="-arch x86_64" \
    ..
make -j${CPUS}
export LIBOMP_AMD64="$BUILD/build_libomp/openmp-${OMPVERSION}.src/build-amd64/runtime/src"

mkdir -p "$BUILD/build_avx2" && cd "$BUILD/build_avx2"
cmake \
    -DCMAKE_INSTALL_PREFIX="${BUILD}/hhsuite" \
    -DCMAKE_BUILD_TYPE=Release \
    -DCHECK_MPI=0 -DHAVE_AVX2=1 \
    -DCMAKE_OSX_ARCHITECTURES=x86_64 \
    -DCMAKE_C_FLAGS="-arch x86_64" -DCMAKE_CXX_FLAGS="-arch x86_64" \
    -DBUILD_SHARED_LIBS=OFF -DCMAKE_FIND_LIBRARY_SUFFIXES=".a" \
    -DOpenMP_C_FLAGS="-Xpreprocessor -fopenmp -I${LIBOMP_AMD64}" -DOpenMP_C_LIB_NAMES=omp -DOpenMP_CXX_FLAGS="-Xpreprocessor -fopenmp -I${LIBOMP_AMD64}" -DOpenMP_CXX_LIB_NAMES=omp -DOpenMP_omp_LIBRARY=${LIBOMP_AMD64}/libomp.a \
    "$REPO"
make -j${CPUS}
make install
check_binary "src/hhblits"

export MACOSX_DEPLOYMENT_TARGET=11.0

mkdir -p "$BUILD/build_libomp/openmp-${OMPVERSION}.src/build-arm64" && cd "$BUILD/build_libomp/openmp-${OMPVERSION}.src/build-arm64"
cmake \
    -DLIBOMP_ENABLE_SHARED=OFF \
    -DLIBOMP_INSTALL_ALIASES=OFF \
    -DLIBOMP_ARCH=aarch64 \
    -DCMAKE_OSX_ARCHITECTURES=arm64 \
    -DCMAKE_C_FLAGS="-arch arm64" \
    -DCMAKE_CXX_FLAGS="-arch arm64" \
    -DLIBOMP_ASMFLAGS="-arch arm64" \
    ..
make -j${CPUS}
export LIBOMP_AARCH64="$BUILD/build_libomp/openmp-${OMPVERSION}.src/build-arm64/runtime/src"

mkdir -p "$BUILD/build_arm64" && cd "$BUILD/build_arm64"
cmake \
    -DCMAKE_BUILD_TYPE=Release \
    -DCHECK_MPI=0 -DHAVE_ARM8=1 \
    -DCMAKE_OSX_ARCHITECTURES=arm64 \
    -DCMAKE_C_FLAGS="-arch arm64" -DCMAKE_CXX_FLAGS="-arch arm64" \
    -DBUILD_SHARED_LIBS=OFF -DCMAKE_FIND_LIBRARY_SUFFIXES=".a" \
    -DOpenMP_C_FLAGS="-Xpreprocessor -fopenmp -I${LIBOMP_AARCH64}" -DOpenMP_C_LIB_NAMES=omp -DOpenMP_CXX_FLAGS="-Xpreprocessor -fopenmp -I${LIBOMP_AARCH64}" -DOpenMP_CXX_LIB_NAMES=omp -DOpenMP_omp_LIBRARY=${LIBOMP_AARCH64}/libomp.a \
    "$REPO"
make -j${CPUS}
check_binary "src/hhblits"

for BINARY in "${BUILD}/hhsuite/bin/"*; do
    BINARY_NAME="$(basename "${BINARY}")"
    case "${BINARY_NAME}" in
        *ffindex*)
          lipo -create \
              -arch x86_64 "$BUILD/build_avx2/lib/ffindex/src/${BINARY_NAME}" \
              -arch arm64 "$BUILD/build_arm64/lib/ffindex/src/${BINARY_NAME}" \
              -output "$BUILD/hhsuite/bin/${BINARY_NAME}"
          ;;
        *)
          lipo -create \
              -arch x86_64 "$BUILD/build_avx2/src/${BINARY_NAME}" \
              -arch arm64 "$BUILD/build_arm64/src/${BINARY_NAME}" \
              -output "$BUILD/hhsuite/bin/${BINARY_NAME}"
          ;;
    esac
done

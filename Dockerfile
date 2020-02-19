ARG NAMESPACE=
FROM debian:stable-slim as qemu-downloader
ARG NAMESPACE
RUN if [ X"$NAMESPACE" != X"" ]; then \
        apt-get update && apt-get install -y wget && rm -rf /var/lib/apt/lists/*; \
    fi; \
    if [ X"$NAMESPACE" = X"ppc64le/" ]; then \
        wget -nv -O /usr/bin/qemu-ppc64le-static https://github.com/multiarch/qemu-user-static/releases/download/v4.2.0-4/qemu-ppc64le-static; \
        chmod +x /usr/bin/qemu-ppc64le-static; \
    fi; \
    if [ X"$NAMESPACE" = X"aarch64/" ]; then \
        wget -nv -O /usr/bin/qemu-aarch64-static https://github.com/multiarch/qemu-user-static/releases/download/v4.2.0-4/qemu-aarch64-static; \
        chmod +x /usr/bin/qemu-aarch64-static; \
    fi; \
    touch /usr/bin/dummy_copy

FROM ${NAMESPACE}debian:stable-slim as builder
ARG NAMESPACE
COPY --from=qemu-downloader /usr/bin/dummy_copy /usr/bin/qemu-ppc64le-static* /usr/bin/
COPY --from=qemu-downloader /usr/bin/dummy_copy /usr/bin/qemu-aarch64-static* /usr/bin/

RUN apt-get update && apt-get install -y \
    build-essential cmake xxd ninja-build \
 && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/hh-suite
ADD . .

WORKDIR /opt/hh-suite/build
RUN cmake -DHAVE_SSE2=1 -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/usr/local/hh-suite ..
RUN make -j $(nproc --all) install

FROM ${NAMESPACE}debian:stable-slim
ARG NAMESPACE
MAINTAINER Milot Mirdita <milot@mirdita.de>
COPY --from=qemu-downloader /usr/bin/dummy_copy /usr/bin/qemu-ppc64le-static* /usr/bin/
COPY --from=qemu-downloader /usr/bin/dummy_copy /usr/bin/qemu-aarch64-static* /usr/bin/

RUN apt-get update && apt-get install -y \
    libstdc++6 libgomp1 \
 && rm -rf /var/lib/apt/lists/*

COPY --from=builder /usr/local/hh-suite /usr/local/hh-suite

ENV HHLIB=/usr/local/hh-suite
ENV PATH="/usr/local/hh-suite/bin:/usr/local/hh-suite/scripts:${PATH}"

CMD ["hhblits"]

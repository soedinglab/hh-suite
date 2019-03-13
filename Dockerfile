FROM debian:stable-slim as builder

RUN apt-get update && apt-get install -y \
    build-essential cmake xxd ninja-build \
 && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/hh-suite
ADD . .

WORKDIR /opt/hh-suite/build
RUN cmake -G Ninja -DHAVE_SSE2=1 -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/usr/local/hh-suite ..
RUN ninja && ninja install

FROM debian:stable-slim
MAINTAINER Milot Mirdita <milot@mirdita.de>

RUN apt-get update && apt-get install -y \
    libstdc++6 libgomp1 \
 && rm -rf /var/lib/apt/lists/*

COPY --from=builder /usr/local/hh-suite /usr/local/hh-suite

ENV HHLIB=/usr/local/hh-suite
ENV PATH="/usr/local/hh-suite/bin:/usr/local/hh-suite/scripts:${PATH}"

CMD ["hhblits"]


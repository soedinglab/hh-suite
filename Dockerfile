FROM alpine:latest as builder

RUN apk add --no-cache gcc g++ cmake musl-dev vim ninja

WORKDIR /opt/hh-suite
ADD . .

WORKDIR /opt/hh-suite/build
RUN cmake -G Ninja -DHAVE_SSE2=1 -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/usr/local/hh-suite ..
RUN ninja && ninja install

FROM alpine:latest
MAINTAINER Milot Mirdita <milot@mirdita.de>
RUN apk add --no-cache libstdc++ libgomp

COPY --from=builder /usr/local/hh-suite /usr/local/hh-suite

ENV HHLIB=/usr/local/hh-suite
ENV PATH="/usr/local/hh-suite/bin:/usr/local/hh-suite/scripts:${PATH}"

CMD ["hhblits"]

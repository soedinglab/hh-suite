FROM alpine:latest as builder

RUN apk add --no-cache gcc g++ cmake musl-dev vim git ninja

WORKDIR /opt/hh-suite
ADD . .

WORKDIR /opt/hh-suite/build
RUN cmake -G Ninja -DHAVE_SSE2=1 -DHAVE_MPI=0 -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/opt/hh-suite ..
RUN ninja && ninja install

WORKDIR /opt/hh-suite
RUN rm -rf /opt/hh-suite/build

FROM alpine:latest
MAINTAINER Milot Mirdita <milot@mirdita.de>
RUN apk add --no-cache bash grep libstdc++ libgomp

COPY --from=builder /opt/hh-suite /opt/hh-suite

ENV HHLIB=/opt/hh-suite
ENV PATH="/opt/hh-suite/bin:/opt/hh-suite/scripts:${PATH}"

CMD ["hhblits"]

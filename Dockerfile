FROM ubuntu:24.04 AS runtime

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        bash \
        build-essential \
        ca-certificates \
        cmake \
        libboost-dev \
        libboost-program-options-dev \
        libeigen3-dev \
        libtbb-dev \
        ninja-build \
        pkg-config \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/eulerian-SPH

ENV ROOT_DIR=/opt/eulerian-SPH \
    BUILD_DIR=/opt/eulerian-SPH/build \
    BUILD_TYPE=Release \
    CMAKE_GENERATOR=Ninja \
    CASE_NAME=2d_eulerian_taylor_green_LG

COPY . /opt/eulerian-SPH
COPY docker/entrypoint.sh /usr/local/bin/eulerian-sph-entrypoint

RUN sed -i 's/\r$//' /usr/local/bin/eulerian-sph-entrypoint \
    && chmod +x /usr/local/bin/eulerian-sph-entrypoint

WORKDIR /opt/eulerian-SPH

ENTRYPOINT ["eulerian-sph-entrypoint"]
CMD []

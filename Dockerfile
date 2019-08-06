#chrisw 20190709
#build an image for running LURE

FROM ubuntu:18.04

ENV tbb_os linux

RUN ["mkdir", "-p", "/dockerfile_scripts/"]
COPY ./dockerfile_scripts /dockerfile_scripts

RUN ["ln", "-s", "/data/repos", "/root/repos"]

RUN ["bash", "/dockerfile_scripts/install_apt_stuff.sh"]
RUN ["bash", "/dockerfile_scripts/install_r_stuff.sh"]

RUN ["mkdir", "-p", "/lure_scripts/"]
COPY ./scripts /lure_scripts

COPY ./running_LURE_from_TACKLE_BOX.md /

ENTRYPOINT ["cat", "/running_LURE_from_TACKLE_BOX.md"]

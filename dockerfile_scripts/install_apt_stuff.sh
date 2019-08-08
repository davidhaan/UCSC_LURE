#!/bin/bash

#chrisw 20190709
#install various packages for running LURE

apt update
apt upgrade -y
apt dist-upgrade -y
apt autoremove -y


# to avoid input prompt for Configuring keyboard-configuration
DEBIAN_FRONTEND=noninteractive apt install keyboard-configuration


# installs add-apt-repository
apt install -y software-properties-common

# x11 stuff
apt install -y xorg xvfb xauth xfonts-base

# install INTEL MKL library
git clone https://github.com/eddelbuettel/mkl4deb.git
bash ./mkl4deb/script.sh

# apt install -y libcurl4-openssl-dev libssl-dev libxml2-dev

# required for Similarity R package
apt install -y libtbb-dev
# export tbb_os=linux


apt install -y python3-pip
pip3 install matplotlib numpy pandas argparse


apt install -y default-jre


# install R3.6, which is not yet in default apt repo.
# https://stackoverflow.com/questions/55929757/installing-r-3-6-on-ubuntu-disco-19-04
# to avoid input prompt for package tzdata
#DEBIAN_FRONTEND=noninteractive apt install -y r-base
add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/'
apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
apt update
DEBIAN_FRONTEND=noninteractive apt install -y r-base-dev



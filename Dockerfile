FROM continuumio/miniconda3

WORKDIR /usr/src/app
COPY ./ ./

RUN mkdir $HOME/mosek
COPY ./mosek.lic /root/mosek/mosek.lic

# Xvfb and other utils
RUN apt-get install -yq --no-install-recommends \
    xvfb \
    x11-utils \
    libx11-dev \
    qt5-default \
	nano \
	vim \
	blender \
    && apt-get clean

ENV DISPLAY=:99

RUN conda env update --name base --file environment.yml

EXPOSE 8080
# The code to run when container is started:
#ENTRYPOINT ["/bin/bash"]
FROM python:slim-stretch
WORKDIR /home
RUN apt-get update -y
RUN DEBIAN_FRONTEND=noninteractive apt-get -yq install libgfortran3
RUN pip3 install sklearn pandas seaborn
ADD ./XRRpred /home/XRRpred
ENV PATH "/home/XRRpred/asaquick/bin:$PATH"
ENTRYPOINT  ["/bin/bash"]




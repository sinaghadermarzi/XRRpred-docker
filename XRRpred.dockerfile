FROM python:slim-stretch
WORKDIR /home
RUN apt-get update -y
RUN DEBIAN_FRONTEND=noninteractive apt-get -yq install libgfortran3
RUN pip3 install --upgrade pip
RUN pip3 install scikit-learn pandas seaborn
ADD ./XRRpred /home/XRRpred
RUN chmod -R 777 /home/XRRpred
ENV PATH "/home/XRRpred/asaquick/bin:$PATH"
ENTRYPOINT  ["/bin/bash"]




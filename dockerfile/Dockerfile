# To run, put this in its own folder, run on terminal:
# docker build --tag idefined .
# This builds the docker image file. You need internet connection for
# downloading of files.
# Next, run:
# docker run -it --publish 8000:8000 idefined
# This should create a docker container and send you into the /home/app/src dir.
# Next, from the container, run:
# /opt/conda/bin/python manage.py runserver 0.0.0.0:8000
# Navigate to your browser, and type localhost:8000

FROM ubuntu
RUN apt update
RUN apt -y upgrade
RUN apt -y install curl gpg
RUN curl https://repo.anaconda.com/pkgs/misc/gpgkeys/anaconda.asc | gpg --dearmor > conda.gpg
RUN install -o root -g root -m 644 conda.gpg /usr/share/keyrings/conda-archive-keyring.gpg
RUN echo "deb [arch=amd64 signed-by=/usr/share/keyrings/conda-archive-keyring.gpg] https://repo.anaconda.com/pkgs/misc/debrepo/conda stable main" > /etc/apt/sources.list.d/conda.list
RUN apt update
RUN apt -y install gcc git
ENV TZ=Asia/Singapore
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone
RUN apt -y install libmpich-dev
RUN curl -s https://packagecloud.io/install/repositories/github/git-lfs/script.deb.sh | bash
RUN apt-get install git-lfs
RUN apt install conda
SHELL ["/opt/conda/bin/conda", "run", "-n", "base", "/bin/bash", "-c"]
RUN conda install python
RUN echo Y | pip install --default-timeout=100 django scikit-learn numpy scipy bokeh matplotlib pandas boto3
RUN export python=/opt/conda/bin/python
ARG CACHEBUST=1
RUN echo "$CACHEBUST"

RUN git clone https://github.com/MelvinYin/Descriptor_Generator.git --depth 1 /home/app
WORKDIR /home/app/src/preprocess/search/source
RUN gcc -o search searchPSSM.c -lm
WORKDIR /home/app/src/preprocess/converge
RUN mpicxx -std=c++11 -Iinclude -o calculator main.cpp
WORKDIR /home/app
RUN mv /home/app/src/preprocess/search/source/search /home/app/src/preprocess/search/search


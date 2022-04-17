
FROM continuumio/miniconda3
WORKDIR /vendorbot_container

COPY env.yml .

RUN conda update -n base -c defaults conda

# install main bot dependencies
RUN conda env create -f env.yml

# manually install node 15 for dev mode nodemon to work
RUN apt-get update -y \
  && apt-get install -y curl \
  && curl -sL https://deb.nodesource.com/setup_16.x | bash - \
  && apt-get install -y nodejs \
  && curl -L https://www.npmjs.com/install.sh | sh

COPY package*.json ./
RUN npm install

# manually install mongodb to access mongodump: https://developer.mongodb.com/community/forums/t/mongodb-org-server-postinst-systemctl-not-found/9232/7
RUN apt-get update -y \
  && apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv 4B7C549A058F8B6B \
  && echo "deb [arch=amd64] http://repo.mongodb.org/apt/ubuntu bionic/mongodb-org/4.2 multiverse" | tee /etc/apt/sources.list.d/mongodb-org-4.2.list \
  && apt-get update -y
RUN ln -T /bin/true /usr/bin/systemctl && apt-get update && apt-get install -y mongodb-org && rm /usr/bin/systemctl

RUN apt-get install -y zip unzip


# RUN conda update -n base -c defaults conda
# RUN conda install -c anaconda capturer

#clone and local install mongo-rdkit repository
RUN git clone https://github.com/rdkit/mongo-rdkit /mongo-rdkit \
  &&  chmod -R 777 /mongo-rdkit && cd /mongo-rdkit && conda run -n vendorbot_env pip install -e .

# copy all other files for bot
COPY . .

# ENV PATH /opt/conda/envs/vendorbot_env/bin:$PATH
# RUN /bin/bash -c "source activate vendorbot_env"

# Make RUN commands use the new environment: -- activate env doesnt work in miniconda!!!
# The SHELL instruction allows the default shell used for the shell form of commands to be overridden. The default shell on Linux is ["/bin/sh", "-c"], and on Windows is ["cmd", "/S", "/C"]. The SHELL instruction must be written in JSON form in a Dockerfile.
# RUN chown -R user /opt/conda/envs
# RUN chmod -R 777 /opt/conda/envs


# Make RUN commands use the new environment: -- activate env doesnt work in miniconda!!!
# The SHELL instruction allows the default shell used for the shell form of commands to be overridden. The default shell on Linux is ["/bin/sh", "-c"], and on Windows is ["cmd", "/S", "/C"]. The SHELL instruction must be written in JSON form in a Dockerfile.
# RUN chown -R user /opt/conda/envs
# RUN chmod -R 777 /opt/conda/envs
#RUN chown -R $USER:$USER /opt/conda/envs

SHELL ["conda", "run", "-n", "vendorbot_env", "/bin/bash", "-c"] 

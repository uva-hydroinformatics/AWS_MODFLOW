FROM ubuntu

WORKDIR /app

ADD . /app


RUN apt-get update

RUN apt-get build-dep -y python-gdal
RUN apt-get install -y python python-dev python-distribute python-numpy gdal-bin libgdal-dev python-pip gcc g++ python-rasterio python-gdal python3-gdal 
RUN pip install flopy
RUN pip install fiona
RUN apt-get upgrade -y

#RUN pwd

#ENTRYPOINT ["ls", "-la", "-R", "/"]

#RUN python build_modflow.py
CMD [ "python", "./build_modflow.py" ]







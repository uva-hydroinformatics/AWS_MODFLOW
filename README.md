# AWS MODFLOW

The code uses Python Flask to trigger the processing of input files from Hydroshare. By sending an HTTP GET request with the Hydroshare ID of the input files it will download them and run the engine.   

To setup this Python Flask webserver on an Ubuntu System follow the steps below:

First obtain the scripts from this repository:  

```git clone https://github.com/uva-hydroinformatics-lab/AWS_MODFLOW.git```  

 On a fresh ubuntu instance install nginx, python, and gunicorn:  
 ```shell
 sudo apt-get install -y python python-pip nginx gunicorn```

 Install the required python packages:
 ```
 pip install flask hs_restclient numpy fiona rasterio flopy```  

 Setup Nginx:
 ```shell
sudo /etc/init.d/nginx start
sudo rm /etc/nginx/sites-enabled/default
sudo touch /etc/nginx/sites-available/flask_project
sudo ln -s /etc/nginx/sites-available/flask_project /etc/nginx/sites-enabled/flask_project```  

 Then edit the config file
 ``` shell
 sudo vim /etc/nginx/sites-enabled/flask_project```  

 ```
 server {
    location / {
        proxy_pass http://localhost:8000;
        proxy_set_header Host $host;
        proxy_set_header X-Real-IP $remote_addr;
    }
}```  

Restart Nginx:
``` shell
sudo /etc/init.d/nginx restart
```  

Start the scrript using gunicorn:
```
cd AWS_MODFLOW
gunicorn app:app -b localhost:8000
```  

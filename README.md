# DEFINED-Proteins
<p>
Defined proteins is now live, at http://definedproteins.com

Note the http and not https, you may have to manually type it in.

Otherwise, either:
1. docker pull melvinyin/definedproteins  
OR 1. Copy the dockerfile/Dockerfile here into a local .txt file named
 Dockerfile, and docker build -t definedproteins --build-arg CACHEBUST=$(date +%s ) .  
Remember the dot at the end.   
2. Because the website loads by default into port 80, you need to kill any
 processes using that port first. sudo lsof -i :80 => sudo kill YOUR_PID   
3. If lsof gives nothing, ignore the kill step. This is usually your default
nginx or apache server intro page, so it's fine. If you're running something
 useful there you know not to do this.[1]  
4. docker run --publish 80:80 melvinyin/definedproteins /bin/bash ./docker_build.sh  
4. If docker gives permission denied, instead of doing sudo docker, grant
 permission to docker instead, via sudo groupadd docker => sudo usermod -aG
  docker $USER => sudo reboot [THIS WILL REBOOT YOUR COMPUTER]. sudo docker
   gives weird errors sometimes. 
5. Navigate to http://localhost in your browser, it should work. Ignore
 security warning.   
 
To run in background in 4, add -d.

If you don't want to use docker, see dockerfile for build instructions. A
 possible set of instructions are:  

sudo apt update  
sudo apt upgrade  
sudo apt install gcc libmpich-dev git    
curl -s https://packagecloud.io/install/repositories/github/git-lfs/script.deb.sh | sudo bash    
sudo apt-get install git-lfs    
For python, either have your own venv activated or, to use conda:  
sudo apt install curl gpg  
curl https://repo.anaconda.com/pkgs/misc/gpgkeys/anaconda.asc | gpg --dearmor > conda.gpg  
install -o root -g root -m 644 conda.gpg /usr/share/keyrings/conda-archive-keyring.gpg  
echo "deb [arch=amd64 signed-by=/usr/share/keyrings/conda-archive-keyring.gpg] https://repo.anaconda.com/pkgs/misc/debrepo/conda stable main" > /etc/apt/sources.list.d/conda.list  
sudo apt update  
sudo apt install conda  
Here, either activate your conda env via:
/opt/conda/bin/conda init bash OR /your/path/to/bin/conda init bash  
. ~/.bashrc  
conda activate base  
conda install python  
echo Y | pip install django scikit-learn numpy scipy bokeh matplotlib pandas
 boto3  
If you're running your own python venv, only run the last line.  
Next, you might have to build the search and converge binaries. After you
 have git cloned:
(I know it's bad to have to cd and therefore set PWD to there, yes it
 implicitly requires PWD to be there, this is old code)  
cd /your/path/to/Defined_Proteins/src/preprocess/converge
mpicxx -std=c++11 -Iinclude -o calculator main.cpp  
cd /your/path/to/Defined_Proteins/src/preprocess/search/source  
gcc -o search searchPSSM.c -lm  
cd ..  
mv ./source/search ./search  
cd ../../..  
(You should be in root now, where src is a directory same level as you and
 you're not in src)  
python download_data_folder.py  

Finally, to launch the program:  
python manage.py runserver 0.0.0.0:8000

Navigate to your browser's localhost:8000 and it should work. 

[1] portno is hard-coded in ./docker_build.sh, you probably need to manually
 docker run -it and change that, sorry. Remember to docker commit/save and
  use that image, and to change port mapping in --publish, in run step.
   
</p>

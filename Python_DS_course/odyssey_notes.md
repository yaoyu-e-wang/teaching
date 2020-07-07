### This guide describes how to setup the proper computing environment for the Introduction to Python for Genomic Analysis course

The guide is based on the general instructions here https://docs.rc.fas.harvard.edu/kb/jupyter-notebook-server-on-odyssey/

---

Assuming you have an account on Odyssey with Research Computing, you should be able to SSH into the server using your login credentials.

From there, log onto a compute node:
```
srun --pty -p shared -t 0-08:00 --mem 6000 /bin/bash
```
Here we use the interactive login flag (`--pty`) on the `shared` queue, ask for 8 hours of time, and request 6GB of RAM.  4-8GB should be sufficient for the tasks in this course.

After logging in, run
```
module load Anaconda3/5.0.1-fasrc02
```
This loads Anaconda, which is a popular package manager.  Anaconda handles the various dependencies between software and installs the software you request.

To build a scientific Python environment that we can reuse, we will create an anaconda "environment": 
```
conda create -y -n jupyter_env python=3.6
```
Typically, Anaconda will ask you before installing anything; the `-y` tells Anaconda that the answer to any questions is "yes".  We named the environment as `jupyter_env`, but feel free to change that to anything you want.  Finally, we tell Anaconda that we want to use Python 3.6 for this environment.

Next, we "activate" the environment and install the required packages (substitute the name for your environment, instead of `jupyter_env`, as necessary):
```
source activate jupyter_env
conda install -y numpy pandas matplotlib scipy scikit-learn jupyterlab nodejs ipympl seaborn
```
(This will typically take a bit of time to complete).  Finally, to enable some interactive plot features, we install an "extension" for Jupyter: 
```
jupyter labextension install @jupyter-widgets/jupyterlab-manager jupyter-matplotlib
```
This will also take some time...

---

### Starting the notebook

To use the notebook, we have to effectively expose a port which will allow you to communicate with Jupyter through the web browser on your local machine.

First, we need to determine a port that is "open" and available to use.  To do that, run:
```
for myport in {6818..11845}; do ! nc -z localhost ${myport} && break; done
echo $myport
```
(The second command should print something like 6819 (or some similar number).  That is the port that will eventually be exposed.

Next, we run
```
echo "ssh -NL $myport:$(hostname):$myport $USER@login.rc.fas.harvard.edu"
```
which will print something like (for me, given that Odyssey placed me on the holy7c18412 node):
```
ssh -NL 6819:holy7c18412.rc.fas.harvard.edu:6819 blawney@login.rc.fas.harvard.edu
```
**We will be using this line of text in a moment, but not yet**.

**Start jupyter with**
```
jupyter lab --ip='0.0.0.0' --no-browser --port=$myport
```
Note that we are using Jupyter "lab", which is practically the same as jupyter notebook, but with a slightly more updated/friendly user-interface.

After a moment, it should print some text similar to
```
[I 12:16:54.332 LabApp] JupyterLab extension loaded from /n/home05/blawney/.conda/envs/test_jupyter_build/lib/python3.6/site-packages/jupyterlab
[I 12:16:54.332 LabApp] JupyterLab application directory is /n/home05/blawney/.conda/envs/test_jupyter_build/share/jupyter/lab
[I 12:16:54.340 LabApp] Serving notebooks from local directory: /n/home05/blawney
[I 12:16:54.340 LabApp] The Jupyter Notebook is running at:
[I 12:16:54.340 LabApp] http://holy7c18412.rc.fas.harvard.edu:6819/?token=8f4acf37ed4ec7fa2f420d8074d597dc977d979d4c374b8a
[I 12:16:54.340 LabApp]  or http://127.0.0.1:6819/?token=8f4acf37ed4ec7fa2f420d8074d597dc977d979d4c374b8a
[I 12:16:54.340 LabApp] Use Control-C to stop this server and shut down all kernels (twice to skip confirmation).
[C 12:16:54.386 LabApp] 
    
    To access the notebook, open this file in a browser:
        file:///n/home05/blawney/.local/share/jupyter/runtime/nbserver-331980-open.html
    Or copy and paste one of these URLs:
        http://holy7c18412.rc.fas.harvard.edu:6819/?token=8f4acf37ed4ec7fa2f420d8074d597dc977d979d4c374b8a
     or http://127.0.0.1:6819/?token=8f4acf37ed4ec7fa2f420d8074d597dc977d979d4c374b8a
```

At this point, Jupyter is running, but we cannot access it yet.  To do that, open *another* terminal window and copy/paste that SSH command in (the one that was like `ssh -NL 6819:holy7c18412.rc.fas.harvard.edu:6819 blawney@login.rc.fas.harvard.edu`).  It will ask for your username and password.  If your login succeeds, nothing will happen, almost like it is not responding.  That is expected.  This command establishes the "bridge" that connects the Jupyter notebook on Odyssey to your local machine.



Now open your web browser and type the address that appeared when Jupyter was started (e.g. `http://127.0.0.1:6819/?token=8f4acf37ed4ec7fa2f420d8074d597dc977d979d4c374b8a` )

This *should* open the Jupyter lab in your browser.

---

### Testing the installation

Here we test that the interactive plotting features installed correctly.  In the Jupyter lab browser, open a new "notebook" from the Launcher.  Then copy/paste the following into the notebook: 

```
# Load the libraries
%matplotlib ipympl
import matplotlib.pyplot as plt
import numpy as np

# number of points:
n = 25

# centers of point clouds:
xc1 = 10
yc1 = 10
xc2 = 15
yc2 = 15

# constants for defining a plane:
a = 1.5
b = 2.0
c = 10

# simulate random points from two sources:
x1 = np.random.normal(loc=xc1, size=n)
y1 = np.random.normal(loc=yc1, size=n)
z1 = c+a*x1+b*y1 + np.random.normal(loc=1, size=n)

x2 = np.random.normal(loc=xc2, size=n)
y2 = np.random.normal(loc=yc2, size=n)
z2 = c+a*x2+b*y2 + np.random.normal(loc=1, size=n)

# points for interpolating surface plot (the plane)
xx, yy = np.meshgrid(np.arange(5,20), np.arange(5,20))
zz = c+a*xx+b*yy

# make the plot
fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(111, projection='3d')
ax.scatter(x1,y1,z1)
ax.scatter(x2,y2,z2)
ax.plot_surface(xx,yy,zz, alpha=0.5)
```

Hit `Shift+Enter` to run the cell.  You should see a 3-dimensional plot with a plane and some random points.  Using your mouse, you should be able to drag the plot to rotate the view.  

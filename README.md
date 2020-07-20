# DCJUC

# 1. Compile and Run DCJUC V1
## 1.1 Compile
before compile, gsl must be installed
how to compile:
```
make
```
## 1.2 Run
how to run (options are in the main.c code):
```
cd main
./DCJUC
```

**IMPORTANT: V1 has been deprecated, it's only for reference, please use V2 which is in folder DCJUC_V2**

# 2. Compile and Run DCJUC V2
## 2.1 Compile 

### 2.1.1 Prerequisite

**MAC OS**

Compilers and configuration tools (you might need to edit configure.ac):

```
brew install gcc@9
brew install autoconf automake libtool
ln -s /usr/local/bin/glibtoolize /usr/local/bin/libtoolize
```

**Ubuntu Linux**

```
sudo apt-get install build-essential
sudo apt-get install autoconf automake gdb git libffi-dev zlib1g-dev libssl-dev
sudo apt install libtool
sudo apt-get install gcc-4.8
```


## 2.2 Run

### 2.2.1 Generate data
``
cd shell/insdis
./batch_gen_dist_graph.sh
``

TODO this is a bug
```
vi data/dist/graph/1000_0.1_0.0_0.1_dual_balanced/0.1_0.0_0.1_0 
# append 0 to the first line

```

### 2.2.2 Run Distance

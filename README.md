# sleepy


## Installation

```console
sudo apt update
sudo apt install -y libboost-all-dev

git clone https://github.com/kevinkorfmann/sleepy
cd sleepy
git clone https://github.com/tskit-dev/kastore
git submodule update --init --recursive

mkdir build
cd build
cmake ..
make
```

Include bin directory into your PATH by modiying the .bashrc file
export PATH=$PATH:/home/*username*/projects/sleepy/bin   

Checkout the minimal example in: https://github.com/kevinkorfmann/sleepy-analysis


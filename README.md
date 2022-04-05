# sleepy

Forward weak dormancy simulations with positive selection.   

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

Include bin directory into your PATH by modiying the .bashrc file:
```console
export PATH=$PATH:/home/*username*/projects/sleepy/bin   
```

Checkout the minimal_code_example.ipynb: 
https://github.com/kevinkorfmann/sleepy-analysis/blob/main/minimal_code_example.ipynb

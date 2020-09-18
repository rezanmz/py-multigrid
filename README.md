# About The Project

This project is an easy-to-understand implementation of **R. Beck's Algebraic Multigrid (AMG) method** published in _December 1999_.

I tried to include as much comment as possible in the code to make it readable.

_Note: This code is written in pure python and written with the intention of understanding and researching AMG methods. If you need a more robust and a more efficient implementation, check out [PyAMG](https://github.com/pyamg/pyamg)._

You can find the original paper at [https://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.569.7668&rep=rep1&type=pdf](https://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.569.7668&rep=rep1&type=pdf)

# Getting Started

To get a local copy of the project, clone this repository on your own computer.

## Installation

1. Clone the repo
```sh
git clone https://github.com/rezanmz/py-multigrid.git
```
2. Change the working directory
```sh
cd py-multigrid
```
3. Create a new virtual environment (Recommended)
```sh
python3 -m venv env && source env/bin/activate
```
4. Install the requirements
```sh
pip install -r requirements.txt
```

# Structure of the project

This project has two main modules:
1. **Numerical Solver** module
2. **Multigrid** module

### Numerical Solver 
This module contains two common and popular iterative solvers, _Jacobi_ and _Gauss Seidel_. These solvers are used as pre-smoother and post-smoother for AMG solver.

### Multgrid
The multigrid module consists of 4 files:
1. **Multigrid**: _The main operations for V-Cycle are handled in this part._
2. **Fine2Coarse**: _This part is responsible for labeling coarse and fine nodes in a given graph._
3. **Prolongator**: _Given coarse and fine nodes found in **Fine2Coarse** part, **Prolongator** finds the prolongation and restriction operator needed to obtain a coarser grid from the fine grid._

# Contributing
Contributions are what make the open source community such an amazing place to be learn, inspire, and create. Any contributions you make are **greatly appreciated**.

1. Fork the Project
2. Create your Feature Branch (`git checkout -b feature/AmazingFeature`)
3. Commit your Changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the Branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

# Contact
If you have any questions of feedbacks, you can get in touch with me:
1. My personal webpage is available at [https://rezanmz.com/](https://rezanmz.com/) and [http://kish.sharif.edu/~reza.namazi/](http://kish.sharif.edu/~reza.namazi/)
2. My email is [rezanmz@ymail.com](mailto:rezanmz@ymail.com)
3. My telegram messenger handle is [@rezanmz](https://t.me/rezanmz)

Project Link: [https://github.com/rezanmz/py-multigrid](https://github.com/rezanmz/py-multigrid)
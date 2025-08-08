# ENGLISH
## Phoenix's Interior Ballistics Solver (PIBS)
A series of Python scripts (with a tkinter Graphical User Interface) for solving the interior ballistics system-of-equations, with provisions for both constrained design and certain optimization calculations.    

## Features
### Calculation from Parameters (Forward Calculation) 
* 0 dimensional interior ballistic SoE in the M.E. Serebryakov formulation.
    
    The interior ballistic problem is formulated in the orthodox, traditional method as is taught in the Soviet Union and in China. The calculation is done in the reduced form for conventional and recoilless guns.

* Efficient Integration

    Integrate up to user specified precision using high-order, adaptive [Runge Kutta Fehlberg](https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta%E2%80%93Fehlberg_method) 7(8)-th order method.

* Supports Shaped Propellant

    Includes sphere, strip, cylinder, and multi-perforated cylinder/prism/rosette.

* Pre-defined Propellant Presets

    Baked in propellant definition for many propellants, sourced from the public literature, with power law, in the form of `a*p^b` burn rate-pressure relation, see `data/propellant.csv`. The burn rate is provided with a fudge factor to allow for granular modification to match experimental result, and can be used to simulate firing at off-normal temperature conditions, if the coefficients are known.


### Constrained Design/Optimization (Backward Calculation):
* Constrained Design conforming to Performance Characteristics

    The required web to achieve a desired peak pressure point can be solved either alone, with option `Lock Gun Length` in the GUI or passing `known_bore=True`, or with desired muzzle velocity, with option `Constrain Design`, in which case the required gun length will be solved simultaneously, for both conventional and recoilless gun.

* Optimization conforming to Performance Characteristics

    From a desired system performance specification, the optimum loading fraction is searched and the one that minimizes total bore volume (the sum of chamber and bore volume) is found via varying the load fraction, with option `Minimize Tube Volume`, for both conventional and recoilless gun.

## How-To
- Use the calculator on Windows:
    - download the latest executable at [latest release](https://github.com/octo-org/octo-repo/releases/latest)
- Development setup:
    - Install Python (>=3.9)
    - Install development dependencies: 
      ```
      pip install .[dev]
      ```
    - (Optional) set up pre-commit hooks
      ```
      exec black . --check
      exec isort . --check
      ```
    - Convert to distributable program via auto-py-to-exe: 
      ```
       auto-py-to-exe -db -c "singlefile.json"
      ``` 

## Contribution
Your contribution is welcomed! Please feel free to raise issues or propose pull requests regarding features or problems encountered. To get started, see the [community documentation](CONTRIBUTING.md)


## Resources Used
* tcl/tk themes used include "awdark" & "awlight" from awthemes: https://wiki.tcl-lang.org/page/awthemes
* monospaced, CJK compatible font: Sarasa-Gothic, https://github.com/be5invis/Sarasa-Gothic


# 中文
## 凤凰内弹道解算器（PIBS）简介
凤凰内弹道计算器是一系列以解火炮，无坐力炮，及高低压发射器的内弹道正反算问题为目的、以Python语言编写并配以基于tk/tcl实现的简易图形界面的计算程序脚本。 

## 功能特性
### 由参数正算
* 实现谢列别梁可夫体系下的0维弹道学
    
    本程序编写紧密参照了我国教科书（如张小兵，金志明《枪炮内弹道学》（2014版，2003版）等），对于常规火炮实现了谢列别梁可夫内弹道体系下，以归一化参数所列内弹道学方程组的解算，并在拉格朗日解的基础上，补充毕杜克（Pidduck）与马蒙陀夫（Mamontov）其他两种0维气动力学梯度，还对于缩膛效应进行了平均补偿与修正。对无座力炮，实现了在金志明拉格朗日解下，以归一化参数所列的无坐力内弹道方程组的解算。

* 高效的积分及求解方法
    
    在参考了前人工作（如王连荣，张佩勤《火炮内弹道计算手册》（1987年）中适用于MZ-80B型微机的BASIC程序，及其它在matlab中编写的内弹道计算程序）的基础上，本程序在编写过程中关注到Python作为解释性语言所固有的计算速度相对较慢的特点，着重优化了内弹道求解的数值算法。对常微分方程的数值求解，采用了自适应的7（8）阶龙格-库塔-费尔伯格（Runge-Kutta-Fehlberg）法。相较于常规的4阶龙格-库塔法，在对方程求解的同时生成残差估计，结合精度需求动态、自控地选取符合要求的步长，克服了固定步长算法中步长选择的盲目性，节约计算步数，实现以较小的计算量快速得出结果。此外，对于一元方程的求零点问题，程序的实现中采用二分与正切混用的德克（Dekker）法。

* 支持不同形状火药
    
    程序内置对于球、柱、管、带、多孔火药的形状函数中，根据几何参数自动计算形状系数、确定相对燃面函数的支持。

* 预设常见火药类型

    程序中收录不少有公开文献的火药燃速，爆温，火药力等数据，并支持用户（在从代码运行时，由于Python打包限制暂不支持可执行文件）补充其他火药，可在`data/propellant.csv`表格中自行添加，程序运行时会按数据定义对应的对象，便于通过图形用户界面直接选择。此外，考虑到内弹道工作实际情况，图形界面另可自定火药燃速修正系数，以便修正批次，药温等因素的影响，使计算结果符合实验。

### 按性能指标反算、优化
* 反算弹道性能对应的设计参数
    
    对于常规火炮及无坐力炮，程序支持按照压强需求反算火药弧厚，以及同时满足弹速要求反算身管长度。

* 优化计算最小膛容解
    
    对于常规火炮以及无座力炮，程序支持计算满足一定压强、弹速性能的最小膛容（含药室容积）解。

## 使用说明
- 下载预先打包为Windows可执行文件的Python代码及运行环境
    - 从[Github发行版直接下载](https://github.com/octo-org/octo-repo/releases/latest)单文件`PIBSva.b.c.d.exe`。
  
    该版本打包Python运行环境及其它依赖库，无法添加火药。单文件为自解压压缩包，每次运行时会解压在系统临时目录，大小约为50MB。

- 在本地Python环境运行代码，或对代码进行再开发
    - 安装[Python运行环境](https://www.python.org/downloads/)。

   
    
 

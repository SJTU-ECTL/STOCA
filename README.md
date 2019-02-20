本程序代码是基于 [[1]](#ref1) 中理论实现的一个样例. 

## (一) 程序的使用

如果直接编译程序，运行的时候请确保目录结构为：

```
./                              # 工作目录
|-- CA-w-h.exe                  # 主程序
|-- mvsis50720.exe              # 主程序会调用它，程序名是 hard code
|-- master.mvsisrc              # mvsis 程序会读取 其实是空文件
|                               # 没有的话会多输出一行字影响主程序读取
|-- pla/                        # 存放结果的目录 目录名是 hard code
|  |-- case-00/                 # 分别存放不同 case 的结果 数字范围是 00 到 99
|  |-- case-01/                 # vector.txt中有多少个case就建多少个文件夹
|  `-- ......                   # 可以使用脚本批量生成
`-- vector.txt                  # 程序的输入 文件名是 hard code
```

### vector.txt 格式

第 1 行: 原始多项式中变量个数 k, 需要的精度 m. 两个整数, 使用空格隔开
第 2 行: 每个变量分别的最高次 d_1, ..., d_k. k 个整数, 用空格隔开
第 3 行: 第一个 case 的 problem vector, d_1 个整数, 用空格隔开
......

从第三行开始就是每个 case 的 problem vector. 

运行程序结束后会在每个 case 的目录下生成结果, 命名为 solution-xxx.pla. 请使用 ABC 进行下一步处理. 

## (二) ABC 处理
ABC 是 Linux 下面的程序, 请在 http://people.eecs.berkeley.edu/~alanmi/abc/ 下载

将 abc.rc, mcnc.genlib, 以及要处理的结果 solution-xxx.pla 和 ABC 程序放在同一个目录下, 运行 ABC 程序. 在 ABC 里输入以下命令

> read_genlib mcnc.genlib
> read solution-xxx.pla
> clp; sop; fx; st; dch; b; map
> write_eqn solution-xxx.eqn

这样就会在目录里生成 solution-xxx.eqn 的结果. 是个文本文档, 记录了得到的逻辑电路. 


## (三) 如何得到 problem vector

### 两种方法: 
#### 第一种方法
使用 MATLAB 程序（见目录 Univariate-Bernstein-Approximation）直接转化为 Bernstein polynomial, 然后对系数进行近似来得到 problem vector.

例: tanh(x)

假设在 MATLAB 程序里选择需要的多项式最高次是 4 (k = 1, d_1 = 4), 运行得到的 Bernstein 系数为

0.0008, 0.2483, 0.5031, 0.6550, 0.7619

于是对应的 Bernstein polynomial 为:

tanh(x) ~= 0.0008*choose(4,0)*x^4 + 0.2483*choose(4,1)*x^3*(1-x) + 0.5031*choose(4,2)*x^2*(1-x)^2 + 0.6550*choose(4,3)*x*(1-x)^3 + 0.7619*choose(4,4)*(1-x)^4

假设需要的精度系数 m = 8, 于是有

0.0008*choose(4,0) ~= 0/256
0.2483*choose(4,1) ~= 254/256
0.5031*choose(4,2) ~= 773/256
0.6550*choose(4,3) ~= 671/256
0.7619*choose(4,4) ~= 195/256

最终得到的 problem vector 就是 V = [0 254 773 671 195]. 如果只有这一个 case, vector.txt 可以写成

1 8
4
0 254 773 671 195

#### 第二种方法
先将函数用其他多项式逼近, 然后再转换为 Bernstein polynomial (其实是转化为一般的 binary combination polynomial, 但 Bernstein polynomial 是等价的, 可以直接使用 [[2]](#ref2) 中的公式)

例: tanh(x)

使用麦克劳林展开式 tanh(x) ~= x - x^3/3 + 2*x^5/15 - 17*x^7/315 + 62*x^9/2835

这里使用一点小变换, 令 t = x^2, f(t) = 1 - t/3 + 2*t^2/15 - 17*t^3/315 + 62*t^4/2835. 那么

tanh(x) = x * f(t)

只需得到 f(t) 的电路, 其输入为 t 的地方与 x 的 squaring unit 相连, 输出用 multiplier 乘以 x.

对于 f(t), 使用 [[2]](#ref2) 中的公式可以得到

b[0] = 1
b[1] = 11/12
b[2] = 77/90
b[3] = 253/315
b[4] = 311/405

那么 f(t) 的 Bernstein polynomial 为

f(t) = b[0]*choose(4,0)*t^4 + b[1]*choose(4,1)*t^3*(1-t) + b[2]*choose(4,2)*t^2*(1-t)^2 + b[3]*choose(4,3)*t*(1-t)^3 + b[4]*choose(4,4)*(1-t)^4

假设 m = 8

b[0]*choose(4,0) ~= 256/256
b[1]*choose(4,1) ~= 939/256
b[2]*choose(4,2) ~= 1314/256
b[3]*choose(4,3) ~= 822/256
b[4]*choose(4,4) ~= 197/256

所以 problem vector 是 V = [256 939 1314 822 197], vector.txt 可以写成

1 8
4
256 939 1314 822 197



## Reference
<a name="ref1">[1]</a> Peng, X., & Qian, W. (2018). Stochastic Circuit Synthesis by Cube Assignment. IEEE Transactions on Computer-Aided Design of Integrated Circuits and Systems, 37(12), 3109-3122.

<a name="ref2">[2]</a> Qian, W., & Riedel, M. D. (2008, June). The synthesis of robust polynomial arithmetic with stochastic logic. In Design Automation Conference, 2008. DAC 2008. 45th ACM/IEEE (pp. 648-653). IEEE.

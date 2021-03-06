#  实验4 数值积分

## 题目1

###  实验要求

![req1](req1.png)

### 算法描述

#### 复合梯形公式求积算法

利用如下公式计算：
$$
\int_{a}^{b} f(x){\rm d}x \approx \frac{h}{2}[f(a) + 2\sum^{n - 1}_{i = 1}f(a + ih) + f(b)]
$$
其中，$h = \frac{b - a}{n}$，$n$为等分份数。

#### 复合辛普森公式求积算法

利用如下公式计算：
$$
\int_{a}^{b} f(x){\rm d}x \approx \frac{h}{6}[f(a) + 4\sum^{n - 1}_{i = 0}f(a + \frac{h}{2} + ih) + 2\sum^{n - 1}_{i = 1}f(a + ih) + f(b)]
$$
其中，$h = \frac{b - a}{n}$，$n$为大区间个数。

#### 龙贝格求积算法

1. 计算$T_0^{(0)} = \frac{b - a}{2}[f(a) + f(b)]$。令$k = 0$。
2. 利用递推公式计算$T_0^{(k)} = \frac{1}{2}T_0^{(k - 1)} + \frac{h}{2}\sum_{i = 0}^{n - 1}f(a + \frac{h}{2} + ih)$。令$k = k + 1$。
3. 利用外推公式计算$T_m^{(m - k)} = \frac{4^mT_{m - 1}^{(m - k + 1)} - T_{m - 1}^{(m - k)}}{4^m - 1}, m = 1, 2, \cdots, k$。
4. 如果$|T_k^{(0)} - T_{k - 1}^{(0)}| \lt \varepsilon$，转5.，否则转2.。
5. $\int_{a}^{b} f(x){\rm d}x \approx T_k^{(0)}$为满足精度要求的积分结果，输出。

实际实现时，每次迭代2. - 4.时，对每个$m$，只需存储$T_m^{(k)}$的最后两个值即可，节约内存空间。

### 程序清单

* `int.cpp`：主要实验代码

### 运行结果

完整输出结果请见`int.txt`。

下述理论积分值为$e - 1 \approx 1.71828182845904509080$。

#### (1) 复合梯形公式以及复合辛普森公式

**复合梯形公式**

余项为$|R_T(f)| = |-\frac{b - a}{12}h^2f''(\eta)|  = \frac{b - a}{12}h^2e^\eta \le \frac{b - a}{12}h^2e^1 = \frac{e}{12n^2}$

令$\frac{e}{12n^2} \le 10^{-6}$，得$n \ge \sqrt{\frac{10^6e}{12}} \approx 475.9$，取$n = 476$

程序输出为

```
  T   = 1.71828246043304599944, diff = 0.00000063197400090864
```

误差为$0.631974 \times 10^{-6}$，与理论值误差不超过$10^{-6}$，理论计算正确。

**复合辛普森公式**

余项为$|R_S(f)| = |-\frac{b - a}{180}(\frac{h}{2})^4f^{(4)}(\eta)|  = \frac{b - a}{180}(\frac{h}{2})^4e^\eta \le \frac{b - a}{180}(\frac{h}{2})^4e^1 = \frac{e}{2880n^4}$

令$\frac{e}{2880n^4} \le 10^{-6}$，得$n \ge (\frac{10^6e}{2880})^{\frac{1}{4}} \approx 5.5$，取$n = 6$，实际上分了$12$个小区间。

程序输出为

```
  S   = 1.71828228843802044423, diff = 0.00000045997897535344
```

误差为$0.459979 \times 10^{-6}$，与理论值误差不超过$10^{-6}$，理论计算正确。

#### (2) 龙贝格求积算法

```
0 1.8591409142295225
1 1.7539310924648253 1.7188611518765928 
2 1.7272219045575166 1.7183188419217472 1.7182826879247577 
3 1.7205185921643018 1.7182841546998968 1.7182818422184403 1.7182818287945305 
  L   = 1.71828182879453050802, diff = 0.00000000033548541722
```

可以看到，$|T_3^{(0)} - T_2^{(0)}| \approx 0.859130 \times 10^{-6} \lt 10^{-6}$时，算法停止。最终结果为$1.71828182879453050802$，与理论值相差$0.000335485 \times 10^{-6}$，远远小于$10^{-6}$。

此时仅分了$2^3 = 8$个等份。

## 题目2

### 实验要求

![req2](req2.png)

### 算法描述

#### 复合高斯公式求积算法

利用如下公式计算：
$$
\int_{a}^{b} f(x){\rm d}x \approx \frac{h}{2}\sum^{n - 1}_{i = 0}[f(a + \frac{h}{2} - \frac{\sqrt{3} h}{6} + ih) + f(a + \frac{h}{2} + \frac{\sqrt{3} h}{6} + ih)]
$$
其中，$h = \frac{b - a}{n}$，$n$为等分份数。

### 程序清单

- `gauss.cpp`：主要实验代码

### 运行结果

本题没有给定精度，这里取为$10^{-6}$。

输出结果请见`gauss.txt`。

$$f^{(4)}(x) = (\frac{4}{1 + x^2})^{(4)} = \frac{1536x^4}{(1 + x^2)^5} - \frac{1152x^2}{(1 + x^2)^4} + \frac{96}{(1 + x^2)^3}$$

$$\max|f^{(4)}(x)| = 96$$

余项为$|R_G(f)| = |\frac{b - a}{4320}h^4f^{(4)}(\eta)| \le \frac{b - a}{4320}h^4 \times 96 =  \frac{1}{45n^4}$

令$\frac{1}{45n^4} \le 10^{-6}$，得$n \ge (\frac{10^6}{45})^{\frac{1}{4}} \approx 12.2$，取$n = 13$。

输出结果为：

```
pi = 3.14159265358979311600
G = 3.14159265368114004602, diff = 0.00000000009134693002
```

与理论值（$\pi$）相差$0.0000913469 \times 10^{-6}$，远远小于$10^{-6}$，这说明上述理论计算给出的界太松了。

经过尝试，取$n = 3$即可让误差小于$10^{-6}$了。

## 体会与展望

本次实验我实现并测试了多种数值积分算法，包括复合梯形公式、复合辛普森公式、龙贝格求积算法以及复合高斯公式。众所周知，求积分自古以来就是一个老大难问题，而数值求积分算法能够用数值的办法，在给定精度要求下，求出一个难以找到原函数的函数的积分值，一般精度还不低。实验中，我体会了龙贝格求积算法收敛速度之快、外推法效果之好。学会数值积分算法，一定能让我今后的学习生活变得更加顺畅。

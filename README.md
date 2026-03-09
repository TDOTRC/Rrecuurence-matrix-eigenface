# Rrecuurence-matrix-eigenface
this repository is created for the paper

现在整体的思路有三种，分别是：

1.计算多种系统的eigenface，并通过这样的eigenface统计系统的特征

2.计算单一recurrence matrix的eigenface，然后根据这样的eigenface 分析背后的特征

3.斜向选取，再计算eigenface 。

针对2.3.两个思路，我认为应该在原有基础上加入以下尝试。

- 确认待测试系统类别

  | 具体类别  | 高斯噪声 | 低频调制 | 漂移调制 |      |
  | :-------- | -------- | -------- | -------- | ---- |
  | 周期函数  |          |          |          |      |
  | 混沌系统1 |          |          |          |      |
  | 混沌系统2 |          |          |          |      |

- 提取其主要特征向量

- 根据主要特征向量进行K聚类算法分类

  

针对1思路，应该改用经过分解后的低维度信息进行K聚类操作

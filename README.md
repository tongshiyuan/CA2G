# CA2G (cDNA position or Amino acid position to genome position)
---

### 简介：
在阅读文献时经常会碰到对一些关键突变的描述不清晰的情况，比如AKT1基因发生了L78T的突变，
有时候甚至不知道转录组的信息，或者cDNA的改变信息，对于少量突变可以采用比如 [Mutalyzer](https://mutalyzer.nl)
等工具进行检索，格式转换，获得genome位置信息。但是一旦数据量太大，手工转换就很麻烦，
开发CA2G，通过cDNA信息或者氨基酸信息快速批量获取可能的genome信息。

### 环境
- 开发环境python 3.9，但是只要是python 3的版本，稍低一些应该也能运行，问题不大。


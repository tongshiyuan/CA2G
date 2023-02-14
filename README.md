# CA2G (cDNA position or Amino acid position to genome position)

## 注意：软件尚未开发完

### 简介：
在阅读文献时经常会碰到对一些关键突变的描述不清晰的情况，比如AKT1基因发生了L78T的突变，
有时候甚至不知道转录组的信息，或者cDNA的改变信息，对于少量突变可以采用比如 [Mutalyzer](https://mutalyzer.nl)
等工具进行检索，格式转换，获得genome位置信息。但是一旦数据量太大，手工转换就很麻烦，
开发CA2G，通过cDNA信息或者氨基酸信息快速批量获取可能的genome信息。

### 环境
- 开发环境python 3.9，但是只要是python 3的版本，稍低一些应该也能运行，问题不大。
- 第三方库：pandas

### 使用
- 快速使用
```shell
python CA2G.py -i example/exp.txt -o example/rst.txt -r hg38
```
- 参数说明
```shell
usage: CA2G.py [-h] -i FILE [-o OUTPUT] [-r REF]

cDNA position or Amino acid position to genome position. v 1.0.0 2023-2-13

optional arguments:
  -h, --help            show this help message and exit
  -i FILE, --file FILE  variants file.
  -o OUTPUT, --output OUTPUT
                        Output of the result. [./<input.file>.pos.txt]
  -r REF, --ref REF     reference version [hg38].
```
1. -f/--file: 输入变异信息文件
2. -o/--output: 输出文件名，默认为当前目录下`<input.file>.pos.txt`，如果文件已存在则会停止运行
3. -r/--ref: 参考基因组，可选为`hg19、hg38`，默认为`hg38`

### 配置文件说明
`DAT`目录下的配置文件来自[ANNOVAR](https://annovar.openbioinformatics.org/en/latest/), `hg19、hg38`对应`-r`参数，
如果有需求，相同格式下可以进行自己配置。

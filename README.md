**dlpipe** 是一个全自动，可断点续传的从 GEO 数据库下载原始 FASTQ 文件的流程。全程仅需提供 GEO 数据库的 GSM 编号，即可全自动搜索下载对应的 FASTQ 文件，支持断点续传、自动重试，适合网络环境一 (la) 般 (ji) 的依赖公开数据苦逼的生信搬运工。

## 依赖安装

推荐使用 [mamba](https://mamba.readthedocs.io/) 安装所需依赖，若你没用安装过 mamba，使用如下命令将其装入 conda 的基础环境

```bash
conda install -n base -c conda-forge mamba
```

使用 mamba 安装所需依赖

```bash
git clone https://github.com/yangmqglobe/dlpipe
cd dlpipe
mamba env create -f environment.yml
conda activate dapipe
```

## 运行流程

流程行运行前，首先需要进行准备工作，创建工作目录

```bash
mkdir -p project/rawdata && cd project/rawdata
```

在目录下创建名为 `gsm.list` 的 GEO 数据库的 GSM 编号列表，每行一个 GSM 编号，形如

```
GSM5104229
GSM5104230
GSM5104231
GSM5104232
GSM5104233
```

此时执行流程

```bash
snakemake --snakefile /path/to/dlpipe/Snakefile --resources network=1 -- dump_all_fastq
```

其中，`--resources network=1` 为设定最大网络并发数，推荐只使用 `1`，更大的并发数并不会显著提高下载速度，仅推荐在某一特定资源速度慢，想优先下载其它资源时使用更大的并发数，但较大的并发数可能导致被封 IP。

其它可选参数有

- `-T` 最大重试次数，如果你的网络非常垃圾，经常下载中断，建议设置为 `0`，即无限次重试。
- `-j` 最大工作线程，可以有效加快 FASTQ 的转化和压缩速度，但是主要限速仍然取决于你的网速，通常不具有很大效果。
- `--config speed_limit=2m` 这个设置项可以对下载进行限速，如果你是在办公室下载，为了保证其他人的上网体验，建议设置，已经贴心为您将该设置定义为仅在早上 9 点到晚上 21 点之间起效，下班后下载也能自动加速。
- `-R metadata` 使用这个配置，会重新整理需要下载的文件的 meta 信息，当你改变了 `gsm.list` 时，会重新评估需要下载的文件列表，建议可以一直使用该选项，即使无任何改动，也不会花费过多的资源。
- 其它 [snakemake 支持的参数请参阅其官方文档](https://snakemake.readthedocs.io/en/stable/executing/cli.html)。

流程将自动在当前工作目录下创建 `fastq/{GSE}/{GSM}/{SRR}` 结构的文件夹，并将对应的 FASTQ 文件放入其中。
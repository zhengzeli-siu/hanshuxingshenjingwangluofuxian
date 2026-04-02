Python/R 代码文件

- readme.txt

- Source_functions_other：包含 FPCA-gY、FPCA-gX、FGLasso 和 PSKL 的 8 个源函数：
  A.prellim.func.R、ADMM.new.R、B.prellim.func.R、FPCA.score.R、ProxAlg_FGM.R、auc.R、base.func.R 和 prec.rec.R。
  所有这些函数均下载自：https://github.com/PercyZhai/FGM_Neighborhood

- Source_functions_npFGM：包含 npFGM 方法的两个源函数：
  main1.R 和 main2.R

- Simulation：
  文件夹 A、B、C 和 D 分别包含模型 A、B、C 和 D 的仿真代码。
  每个文件夹下又包含三个对应不同维度的子文件夹。
  例如：
  - p30 → p = 30
  - p50 → p = 50
  - p100 → p = 100

- 在模型 A 的子文件夹中（JCGS_Code → Simulation → A → p30），包含 4 个文件：
  classifier.sh、classifier.py、npFGM.R 和 data_generator.R。

- 在模型 B、C 和 D 的每个子文件夹中（例如 JCGS_Code → Simulation → B → p30），包含 5 个文件：
  classifier.sh、classifier.py、npFGM.R、other.R 和 data_generator.R。

- Real_data：包含 5 个文件夹，分别包括真实数据（ADHD）、源文件、R/Python 脚本以及分析结果。
  - data：time.series.ADHD.Rdata 是原始 ADHD 数据集。
  - npFGM：R 脚本 run_ADHD.R 用于执行 npFGM 方法的数据分析。
    Results 子文件夹中包含：
    - ADHD1_npFGM.Rdata：ADHD 组结果
    - ADHD2_npFGM.Rdata：对照组结果
  - fgDNN：包含两个 Python 脚本（fgDNN_ADHD.sh 和 fgDNN_ADHD.py）以及两个子文件夹（get_score 和 results）。
    输出结果保存在 results 文件夹中。
  - Other_alternatives：Results 子文件夹包含 6 个 .Rdata 文件（FPCA-gY、FPCA-gX、FGLasso 和 PSKL 方法的结果）。
    Results Analysis 子文件夹包含 3 个用于分析这些结果的 R 脚本。
    所有这些文件均下载自：https://github.com/PercyZhai/FGM_Neighborhood
  - plot：运行 plot.R 可生成最终图形和分析结果。AAL.info.R 包含图谱（atlas）的坐标信息。

---

如何运行仿真代码：
1. 将 Source_functions_npFGM 和 Source_functions_other 中的所有源文件保存到你的工作目录中。
2. 对于每一种数据生成设定，进入 Simulation 中对应的子文件夹。例如：JCGS_Code → Simulation → A → p30，然后运行 data_generator.R 以生成模型 A 且 p = 30 的数据。
3. 对于 fgDNN 方法：运行 classifier.py（Python 脚本）。bash 脚本 classifier.sh 可用于在 HPC 上提交作业。
4. 对于 npFGM 方法：进入相应子文件夹并运行 npFGM.R。
5. 对于模型 A，替代方法（FPCA-gY、FPCA-gX、FGLasso 和 PSKL）的代码未包含在此处；其结果可从以下地址获取：https://github.com/PercyZhai/FGM_Neighborhood
6. 对于模型 B、C 和 D：运行 other.R 以执行 FPCA-gY、FPCA-gX、FGLasso 和 PSKL 方法。

---

如何运行真实数据分析代码：
1. 进入 Real_data 文件夹。
2. 将 data 中的 time.series.ADHD.Rdata，以及 Source_functions_npFGM 和 Source_functions_other 中的所有源文件保存到你的工作目录中。
3. 对于 npFGM 方法：进入 npFGM 文件夹并运行 run_ADHD.R。
4. 对于 FPCA-gY、FPCA-gX、FGLasso 和 PSKL 方法：结果位于 Other_alternatives → Results 中的 6 个 .Rdata 文件内；分析脚本位于 Other_alternatives → Results Analysis。
5. 对于 fgDNN 方法：
   a. 在 get_score 中运行 ADHD_score.R 以提取分数（保存为 ADHD_1.csv 和 ADHD_2.csv）。
   b. 运行 fgDNN_ADHD.sh，以调用 fgDNN_ADHD.py 并得到结果。
6. 若需生成最终图形：进入 plot 文件夹并运行 plot.R。AAL.info.R 包含图谱坐标信息。

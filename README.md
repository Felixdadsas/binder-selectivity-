# binder-selectivity

> **Do generative protein design models discriminate between homologous targets? BNP vs ANP case study.**

## 1. 问题背景 (Problem)
Brain Natriuretic Peptide (BNP) 与 Atrial Natriuretic Peptide (ANP) 属于高度同源的肽类激素靶标。本库旨在探究主流生成式蛋白质设计模型在处理高度同源靶标时，是否具备真正的选择性区分能力。

## 2. 复现指南 (Reproducibility)
- **输入定义 (`targets/`)**: 标准化的靶标结构 YAML 文件。
- **评测流程 (`notebooks/`)**: 基于 Colab 的端到端推理与打分脚本。

## 3. 免责声明 (Disclaimer)
- **本库内容无实验验证（No experimental validation），所有结论与数据均只限于模型自身的计算打分（In silico scores only）。**

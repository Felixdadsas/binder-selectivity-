# Selectivity Benchmark Notebook Notes

## 运行避坑指南 (Troubleshooting & Pitfalls)
1. **参数设置**: 运行推理时务必加上 `--use_kernels false` 参数，以避免部分云端环境由于缺少特定硬件加速内核导致的编译报错。
2. **环境依赖**: 推荐在 Google Colab 的 GPU 运行实例（如 T4/L4）下执行，确保依赖库版本与输入 YAML 的残基约束对齐。

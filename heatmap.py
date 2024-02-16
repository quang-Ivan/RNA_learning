import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

# 加载数据
data = pd.read_csv('SRR7961226_gene_abundances.tsv', sep='\t')

# 选择要绘制的列，这里我们选择FPKM和TPM作为表达量的代表
# 你可以根据需要选择其中一个，或者如果有多个样本，选择多列
expression_data = data[['Gene ID', 'FPKM', 'TPM']].set_index('Gene ID')
# 对数据根据TPM平均值进行排序
expression_data_sorted = expression_data.sort_values(by='TPM', ascending=False)[0:100]

# 对排序后的数据进行对数转换以改善可视化
log_expression_data_sorted = np.log1p(expression_data_sorted)

# 绘制热图
plt.figure(figsize=(10, 8))
sns.heatmap(log_expression_data_sorted, cmap='viridis')
plt.title('Gene Expression Heatmap (Sorted by TPM)')
plt.ylabel('Gene Name')
plt.xlabel('Expression Measures')
plt.show()

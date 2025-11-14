import pandas as pd

input_file = "/public/compute_dG/result_dG/UniMoMo_pyrosetta_dG.jsonl"
output_file = "/public/compute_dG/result_dG/UniMoMo_sort_pyrosetta_dG.jsonl"  

df = pd.read_csv(input_file, sep="\t")

# 按 foldx_dG 从小到大排序
df_sorted = df.sort_values(by="pyrosetta_dG", ascending=True)

df_sorted.to_csv(output_file, sep="\t", index=False)

print(f" 排序完成！已保存为 {output_file}")

import pandas as pd

# Загрузка данных
data = pd.read_csv("shannon_indices.tsv", sep='\t')

# Вычисление общей статистики
mean_shannon = data['ShannonIndex'].mean()
median_shannon = data['ShannonIndex'].median()
min_shannon = data['ShannonIndex'].min()
max_shannon = data['ShannonIndex'].max()
std_shannon = data['ShannonIndex'].std()

# Вывод общей статистики
print(f"Mean: {mean_shannon:.7f}")
print(f"Median: {median_shannon:.7f}")
print(f"Min: {min_shannon:.7f}")
print(f"Max: {max_shannon:.7f}")
print(f"Std Dev: {std_shannon:.7f}")

# Порог фильтрации (Mean + 3 × Std Dev)
threshold = mean_shannon + 3 * std_shannon

# Фильтрация
filtered_data = data[data['ShannonIndex'] > threshold]
# Печать сколько ортогрупп получилось
print(f"Ортогрупп для анализа: {len(filtered_data)}")

# Запись отфильтрованных данных
filtered_data.to_csv("shannon_indices_filt.tsv", sep='\t', index=False)

print(f"Filtered groups saved to 'shannon_indices_filt.tsv' (threshold: {threshold:.7f})")
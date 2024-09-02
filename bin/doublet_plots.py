#!/usr/bin/env python

import pandas as pd
import matplotlib.pyplot as plt
Data=pd.read_csv('all_doublet_results_combined.tsv',sep='\t')

Data = Data.filter(regex='^(?!.*[s|S]core).*$')
Data = Data.drop(columns=['scrublet__predicted_multiplet'])
# Data= Data.set_index('barcodes')
count_data = Data.apply(pd.Series.value_counts).fillna(0)
count_data = count_data.drop(columns=['barcodes'])

summary_data = count_data.iloc[-2:]  # Assuming the last two rows are the summary
summary_data
# Transpose the data for plotting
summary_data_transposed = summary_data.T

# Plotting
summary_data_transposed.plot(kind='bar', stacked=True)
plt.title('Distribution of Singlets and Doublets Across Different Droplet Types')
plt.xlabel('Droplet Type')
plt.ylabel('Count')
plt.xticks(rotation=45, ha='right')
plt.legend(title='Droplet Type')
plt.tight_layout()
plt.savefig('droplet_type_distribution.png', dpi=300, bbox_inches='tight')

# Show plot
plt.show()
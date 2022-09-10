# Example real data analysis pipeline:

- make sure you have downloaded and replaced the data files in `data` that are stubs
- now run the following commands from the top directory of the github repository:
```
script_analysis/genotype_WGS.sh
script_analysis/err2info.sh
script_analysis/combine_infos.sh
script_analysis/get_metrics.sh
```
- examine the metrics computed in files `data/CAPG*`

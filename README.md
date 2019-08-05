# MassSpecDashBoard

mass spectrometer dashboard for every user to master the status for the big toy. Happy mass spec-ing :)

## How to use:

`python script/massspecdashboard.py your_config_file`

config file format:

```yaml
identifier1: file directory1
identifier2: file directory2
identifier3: file directory3
...
```

for PD_ID, please use result reporter node of PD 2.x. ignore the suffix of PD output. e.g.:

with identifier of `condition 1`:

```
directory/your_filename_Proteins.txt
directory/your_filename_PSM.txt
directory/your_filename_PeptideGroups.txt
```

we should have `condition 1: directory/your_filename` in YAML config file.
